# Exact Diagonalization of the Hubbard Model

I decided to re-explore a [challenging problem](http://qmcchem.ups-tlse.fr/files/caffarel/Hub_Inf_PRL_1994.pdf) first encountered in graduate school as an excuse to experiment with modern web technologies, specifically [WebAssembly](https://webassembly.org/) and [Serverless](https://serverless.com).

Because a large part of the challenge is memory, I started by implementing a version of [this stackoverflow answer](https://stackoverflow.com/a/36345790/8479938) in [AssemblyScript](https://github.com/AssemblyScript/assemblyscript).  I wasn't impressed by the compiler output so I hand tuned the [.wat](http://webassembly.github.io/spec/core/text/index.html) directly because why not.
