(module
(type $tindex (func (param i32) (result i32)))
(type $tmask (func (param i32 i32 i32 i32) (result i32)))
(func $index (export "index") (type $tindex) (param $mask i32) (result i32)
  (local $n i32)
  (local $k i32)
  (local $nmkp1 i32)
  (local $nCk i32)
  (local $nCkm1 i32)
  (local $index i32)
  i32.const 1
  set_local $n        ;; n = 1
  i32.const 1
  set_local $k        ;; k = 1;
  i32.const 1
  set_local $nmkp1    ;; n-k+1 = 1;
  i32.const 0
  set_local $nCk      ;; C(n, k) = 0;
  i32.const 1
  set_local $nCkm1    ;; C(n, k-1) = 1;
  i32.const 0
  set_local $index    ;; index = 0;
  block $break
   loop $repeat
    get_local $mask
    i32.const 1
    i32.and
    if                ;; if (mask & 1)
     get_local $index
     get_local $nCk
     i32.add          ;; index += C(n, k)
    set_local $index
     get_local $mask
     i32.const 1
     i32.eq
     br_if $break     ;; if (mask == 1) break
     get_local $nCkm1
     get_local $n
     i32.mul          ;; C(n, k-1) *= n (and leave on stack)
     get_local $k
     i32.div_u
     set_local $nCkm1 ;; C(n, k-1) *= k
     get_local $nCk
     get_local $n
     i32.mul          ;; C(n, k) *= n (and leave on stack)
     get_local $k
     i32.const 1
     i32.add
     tee_local $k     ;; ++k (and leave on stack)
     i32.div_u
     set_local $nCk   ;; C(n, k) *= k
    else    
     get_local $nCk
     get_local $nCkm1
     i32.add
     set_local $nCk   ;; C(n, k) += C(n, k-1)
     get_local $nCkm1
     get_local $n
     i32.mul          ;; C(n, k-1) *= n (and leave on stack)
     get_local $nmkp1
     i32.div_u
     set_local $nCkm1 ;; C(n, k-1) /= n-k+1
     get_local $nmkp1
     i32.const 1
     i32.add
     set_local $nmkp1 ;; ++(n-k+1)
    end
    get_local $mask
    i32.const 1
    i32.shr_u
    set_local $mask   ;; mask >>= 1
    get_local $n
    i32.const 1
    i32.add
    set_local $n      ;; ++n
    br $repeat
   end
  end
  get_local $index)   ;; return index
(func $mask (export "mask") (type $tmask) (param $index i32) (param $n i32) (param $k i32) (param $nCk i32) (result i32)
  (local $nmk i32)
  (local $mask i32)
  get_local $n
  get_local $k
  i32.sub
  set_local $nmk      ;; n-k
  i32.const 0
  set_local $mask     ;; mask = 0
  block $break
   loop $repeat
    get_local $index
    get_local $nCk
    i32.ge_u
    if                ;; if (index >= C(n, k))
     get_local $mask
     i32.const 1
     i32.or
     set_local $mask  ;; mask |= 1
     get_local $k
     i32.const 1
     i32.eq
     br_if $break     ;; if (k == 1) break
     get_local $index
     get_local $nCk
     i32.sub
     set_local $index ;; index -= C(n, k)
     get_local $nCk
     get_local $k
     i32.mul          ;; C(n, k) *= k (and leave on stack)
     get_local $n
     i32.div_u
     set_local $nCk   ;; C(n, k) /= n
     get_local $k
     i32.const 1
     i32.sub
     set_local $k     ;; --k
    else    
     get_local $nCk
     get_local $nmk
     i32.mul          ;; C(n, k) *= n-k (and leave on stack)
     get_local $n
     i32.div_u
     set_local $nCk   ;; C(n, k) /= n
     get_local $nmk
     i32.const 1
     i32.sub
     set_local $nmk   ;; --(n-k)
    end
    get_local $mask
    i32.const 1
    i32.shl
    set_local $mask   ;; mask <<= 1
    get_local $n
    i32.const 1
    i32.sub
    set_local $n      ;; --n
    br $repeat
   end
  end
  get_local $mask
  get_local $n
  i32.shl)            ;; return mask << n
(func $count (export "count") (type $tindex) (param $mask i32) (result i32)
  get_local $mask
  i32.popcnt))