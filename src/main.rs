use std::collections::TryReserveError;
use std::process::ExitCode;
use std::time::Instant;

use clap::Parser;

use hubbard::basis::Basis;
use hubbard::hamiltonian::Hamiltonian;
use hubbard::lanczos::{ground_state, lanczos, residual};
use hubbard::lattice::{Irrep, Momentum, Symmetry, TiltedSquare};

/// Ground state of the Hubbard model on a tilted square lattice with periodic boundaries.
#[derive(Parser)]
struct Args {
    /// Number of sites, a sum of two squares between 2 and 32
    #[arg(short = 'n', long)]
    sites: usize,
    /// Number of spin up fermions
    #[arg(long)]
    up: u32,
    /// Number of spin down fermions
    #[arg(long)]
    down: u32,
    /// Hopping integral t
    #[arg(short = 't', long, default_value_t = 1.0, allow_negative_numbers = true)]
    hopping: f64,
    /// On-site interaction U
    #[arg(short = 'U', long, allow_negative_numbers = true)]
    interaction: f64,
    /// Stop when the Ritz residual ‖Hψ - Eψ‖ falls below this
    #[arg(long, default_value_t = 1e-8)]
    tolerance: f64,
    #[arg(long, default_value_t = 1000)]
    max_iterations: usize,
    /// Crystal momentum of the symmetry sector
    #[arg(long, value_enum, default_value_t = Momentum::Gamma)]
    momentum: Momentum,
    /// Point group irrep of the symmetry sector
    #[arg(long, value_enum, default_value_t = Irrep::A1)]
    irrep: Irrep,
    /// Try every sector with real characters and report the lowest energy
    #[arg(long)]
    all_sectors: bool,
    /// Ignore symmetries: one sector holding every configuration
    #[arg(long, conflicts_with_all = ["momentum", "irrep", "all_sectors"])]
    full_basis: bool,
    /// Also build the ground state (one more vector) to verify it and measure double occupancy
    #[arg(long)]
    eigenvector: bool,
}

fn main() -> ExitCode {
    match run(Args::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("error: {message}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: Args) -> Result<(), String> {
    let start = Instant::now();
    let lattice = Some(args.sites)
        .filter(|n| (2..=32).contains(n))
        .and_then(TiltedSquare::new)
        .ok_or(format!("{} sites is not a sum of two squares between 2 and 32", args.sites))?;
    if args.up as usize > args.sites || args.down as usize > args.sites {
        return Err(format!(
            "at most {} fermions per spin fit on {} sites",
            args.sites, args.sites
        ));
    }

    println!("lattice    {} sites, tilt {:?}", args.sites, lattice.tilt());

    if args.full_basis {
        let identity = vec![((0..args.sites as u8).collect(), false)];
        solve(&args, &lattice, identity)?.ok_or("no states with these fillings")?;
    } else if args.all_sectors {
        let mut lowest: Option<(f64, String)> = None;
        for momentum in [Momentum::Gamma, Momentum::M] {
            for irrep in [Irrep::A1, Irrep::A2, Irrep::B1, Irrep::B2] {
                if lattice.is_chiral() && matches!(irrep, Irrep::A2 | Irrep::B2) {
                    continue;
                }
                let Ok(group) = lattice.symmetries(momentum, irrep) else { continue };
                println!("sector     {momentum:?} {irrep:?}");
                if let Some(energy) = solve(&args, &lattice, group)?
                    && lowest.as_ref().is_none_or(|(e, _)| energy < *e)
                {
                    lowest = Some((energy, format!("{momentum:?} {irrep:?}")));
                }
            }
        }
        let (energy, sector) = lowest.ok_or("no states with these fillings")?;
        println!("lowest     {energy:.12} in {sector}");
    } else {
        let group = lattice.symmetries(args.momentum, args.irrep)?;
        solve(&args, &lattice, group)?.ok_or("no states of this symmetry with these fillings")?;
    }
    println!("time       {:.2?}", start.elapsed());
    Ok(())
}

/// Ground state energy within the sector spanned by `group`, or `None` if the sector is empty.
fn solve(args: &Args, lattice: &TiltedSquare, group: Vec<Symmetry>) -> Result<Option<f64>, String> {
    let basis = Basis::new(group, args.up, args.down);
    let dimension = basis.dimension();
    let vectors = if args.eigenvector { 3 } else { 2 };
    println!(
        "basis      {dimension} states, {} symmetries, {} for {vectors} vectors",
        basis.group_order(),
        bytes(vectors * dimension * 8)
    );
    if dimension == 0 {
        return Ok(None);
    }

    let h = Hamiltonian::new(basis, &lattice.neighbors(), args.hopping, args.interaction);
    let result = lanczos(&h, args.tolerance, args.max_iterations).map_err(out_of_memory)?;
    let n = args.sites as f64;
    println!("lanczos    {} iterations, residual {:.1e}", result.iterations(), result.residual());
    println!("energy     {:.12} ({:.12} per site)", result.energy, result.energy / n);

    if args.eigenvector {
        let psi = ground_state(&h, &result).map_err(out_of_memory)?;
        let (deviation, energy) = residual(&h, &psi).map_err(out_of_memory)?;
        println!("⟨ψ|H|ψ⟩    {energy:.12}, ‖Hψ - Eψ‖ {deviation:.1e}");
        println!("⟨n↑n↓⟩     {:.12} per site", h.double_occupancy(&psi) / n);
    }
    Ok(Some(result.energy))
}

fn out_of_memory(error: TryReserveError) -> String {
    format!("not enough memory for the Lanczos vectors ({error})")
}

fn bytes(count: usize) -> String {
    let units = ["B", "KiB", "MiB", "GiB", "TiB"];
    let power = ((count.max(1).ilog2() / 10) as usize).min(units.len() - 1);
    format!("{:.1} {}", count as f64 / (1u64 << (10 * power)) as f64, units[power])
}
