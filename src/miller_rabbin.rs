// Power set
// Rabin miller primality test



// Based on miller rabin:
// +1 -> Composite
// -1 -> Probably prime


use rand::Rng;

// Checks if a number is a prime number using miller rabin primality test
pub fn is_prime(n: usize) -> bool {

    let a = rand::thread_rng().gen_range(2..n - 2);

    let mut i = 1_usize;

    loop {
        let mut s = 0;
        let mut d = 0;

        if (n - 1) % (1 << i) == 0{
            s = i;
            d = (n - 1) / (1 << i);

            let mut x = a.pow(d as u32) % n;
            let y = (x * x) % n;

            for i in 0..s {

                if y == 1 && x != 1 && x != n - 1 {
                    return false;
                }

                x = y;
            }

            if y != 1 {
                return false;
            }


            return true;

        }
        
    }
    
}

