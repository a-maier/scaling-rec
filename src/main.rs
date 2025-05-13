use env_logger::Env;
use expression::EvalCountingExpression;
use log::info;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rare::{algebra::{poly::flat::FlatPoly, rat::Rat}, rec::primes::LARGE_PRIMES, traits::{One, TryEval}, Z64};
use rug::Integer;
use seq_macro::seq;

mod expression;
mod parser;

fn main() {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let expr = EvalCountingExpression::<2>::try_from(
        include_str!("../data/coeff_prop_4l")
    ).unwrap();
    info!("Reconstructing coefficient of four-loop propagator");
    cmp_rec(expr);

    let expr = EvalCountingExpression::<2>::try_from(
        include_str!("../data/coeff_prop_4l_mod")
    ).unwrap();
    info!("Reconstructing coefficient of four-loop propagator w/out monomial factors");
    cmp_rec(expr);

    let expr = EvalCountingExpression::<4>::try_from(
        include_str!("../data/aajamp")
    ).unwrap();
    info!("Reconstructing coefficient in two-loop diphoton plus jet amplitude");
    cmp_rec(expr);

    let expr = EvalCountingExpression::<4>::try_from(
        include_str!("../data/aajamp_mod")
    ).unwrap();
    info!("Reconstructing coefficient in two-loop diphoton plus jet amplitude w/out monomial factors");
    cmp_rec(expr);
}

fn cmp_rec<const N: usize>(
    rat: EvalCountingExpression::<N>
) {
    let res_scal = rec_scaling(&rat, Xoshiro256StarStar::seed_from_u64(1));
    rat.print_and_reset_count();
    info!(
        "Reconstructed function has {} independent coefficients",
        res_scal.num().len() + res_scal.den().len() - 1
    );
}

fn rec_scaling<const N: usize>(
    orig: &EvalCountingExpression<N>,
    mut rng: impl Rng,
) -> Rat<FlatPoly<Integer, N>> {
    use rare::rec::rat::thiele_multivar::{Rec, Status};
    use std::ops::ControlFlow::*;

    info!("Starting \"scaling\" reconstruction");
    const P: u64 = LARGE_PRIMES[0];

    let mut rec: Rec<P, N> = Rec::with_random_shift(1, &mut rng);
    let (mut z, mut q_z) = std::iter::repeat_with(|| {
        let z: [Z64<P>; N] = [(); N].map(|_| rng.gen());
        orig.try_eval(&z).map(|q_z| (z, q_z))
    })
        .flatten()
        .next()
        .unwrap();
    while let Continue(Status::Varying(n)) = rec.add_pt(z, q_z).unwrap() {
        z[n] += Z64::one();
        q_z = loop {
            if let Some(val) = orig.try_eval(&z) {
                break val;
            }
            z[n] += Z64::one() ;
        };
    }
    seq!{ M in 0..10 {{
        const P: u64 = LARGE_PRIMES[M];
        loop {
            let z: Z64<P> = rng.gen();
            let z = rec.to_args(z).unwrap();
            let q_z = orig.try_eval(&z).unwrap();

            match rec.add_pt(z, q_z).unwrap() {
                Continue(Status::Scaling) => { },
                Continue(Status::NextMod) => break,
                Continue(Status::Varying(_)) => unreachable!(),
                Break(res) => return res,
            }
        }
    }}}
    panic!("Need more than 10 characteristics!");
}
