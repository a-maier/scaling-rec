use env_logger::Env;
use expression::EvalCountingExpression;
use log::info;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rare::{algebra::{poly::flat::FlatPoly, rat::Rat}, rec::{primes::LARGE_PRIMES, probe::Probe}, traits::{One, TryEval}, Z64};
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
    cmp_rec_all(expr);

    let expr = EvalCountingExpression::<4>::try_from(
        include_str!("../data/aajamp")
    ).unwrap();
    info!("Reconstructing coefficient in two-loop diphoton plus jet amplitude");
    cmp_rec_fast(expr);
}

fn cmp_rec_all<const N: usize>(rat: EvalCountingExpression::<N>) {
    let res_lin = rec_thiele_linear(&rat, Xoshiro256StarStar::seed_from_u64(1));
    rat.print_and_reset_count();

    let res_fast = cmp_rec_fast(rat);
    assert_eq!(res_lin, res_fast);
}

fn cmp_rec_fast<const N: usize>(
    rat: EvalCountingExpression::<N>
) -> Rat<FlatPoly<Integer, N>> {
    let res_hom = rec_homogeneous(&rat, Xoshiro256StarStar::seed_from_u64(1));
    rat.print_and_reset_count();

    let res_scal = rec_scaling(&rat, Xoshiro256StarStar::seed_from_u64(1));
    rat.print_and_reset_count();
    assert_eq!(res_hom, res_scal);
    info!(
        "Reconstructed function has {} independent coefficients",
        res_scal.num().len() + res_scal.den().len() - 1
    );
    res_hom
}

fn rec_thiele_linear<const N: usize>(
    orig: &EvalCountingExpression<N>,
    mut rng: impl Rng,
) -> Rat<FlatPoly<Integer, N>>
{
    use rare::rec::rat::thiele_linear::{Rec, Status};

    info!("Starting Thiele + linear reconstruction");
    const P: u64 = LARGE_PRIMES[0];

    let mut rec = Rec::new(1);
    let (mut z, mut q_z) = std::iter::repeat_with(|| {
        let z: [Z64<P>; N] = [(); N].map(|_| rng.gen());
        orig.try_eval(&z).map(|q_z| (z, q_z))
    })
        .flatten()
        .next()
        .unwrap();
    while let Status::Degree(n) = rec.add_pt(z, q_z) {
        z[n] += Z64::one();
        q_z = loop {
            if let Some(val) = orig.try_eval(&z) {
                break val;
            }
            z[n] += Z64::one() ;
        };
    }
    // TODO: re-use points
    let Status::Rat(mut npts) = rec.status() else {
        panic!("Failed to reconstruct degrees")
    };
    seq!{ M in 0..20 {{
        const P: u64 = LARGE_PRIMES[M];
        let pts = Vec::from_iter(
            std::iter::repeat_with(|| {
                let arg: [Z64<P>; N] = [(); N].map(|_| rng.gen());
                orig.try_eval(&arg).map(|val| Probe{arg, val})
            })
                .flatten()
                .take(npts)
        );
        match rec.add_pts(&pts) {
            Status::Rat(n) => npts = n,
            Status::Done => return rec.into_rat().unwrap(),
            _ => unreachable!("Unexpected reconstruction return value")
        }
    }}}
    let _ = npts;
    panic!("Need more than 20 characteristics!");
}

fn rec_homogeneous<'a, const N: usize>(
    orig: &'a EvalCountingExpression<N>,
    rng: impl Rng,
) -> Rat<FlatPoly<Integer, N>> {
    info!("Starting homogeneous reconstruction");
    match N {
        2 => {
            let rat: &'a EvalCountingExpression<2> = unsafe {
                std::mem::transmute(orig)
            };
            let res = rec_homogeneous_2(rat, rng);
            unsafe{ std::mem::transmute(res) }
        },
        4 => {
            let rat: &'a EvalCountingExpression<4> = unsafe {
                std::mem::transmute(orig)
            };
            let res = rec_homogeneous_4(rat, rng);
            unsafe{ std::mem::transmute(res) }
        },
        _ => unimplemented!("Homogenous reconstruction in {N} variables")
    }
}

fn rec_homogeneous_2(
    orig: &EvalCountingExpression<2>,
    mut rng: impl Rng,
) -> Rat<FlatPoly<Integer, 2>>
{
    use rare::rec::rat::cuyt_lee::{Rec2, Status};
    use rare::rec::rat::finite::cuyt_lee::Needed;
    const NVARS: usize = 2;

    let mut rec = Rec2::with_shift(1, rng.gen());
    seq!{ N in 0..20 {{
        const P: u64 = LARGE_PRIMES[N];
        let z: [Z64<P>; NVARS] = [(); NVARS].map(|_| rng.gen());
        let q_z = orig.try_eval(&z).unwrap();
        rec.add_pt(z, q_z).unwrap();
        loop {
            match rec.status().unwrap() {
                Status::Done => return rec.into_rat().unwrap(),
                Status::NextMod => break,
                Status::Needed(Needed::Pts(pts)) => {
                    let pts = Vec::from_iter(
                        pts.iter()
                            .filter_map(
                                |z| orig.try_eval(z).map(|q_z| (*z, q_z))
                            )
                    );
                    rec.add_pts(&pts).unwrap();
                },
                Status::Needed(Needed::Pt(z)) => {
                    let q_z: Z64<P> = orig.try_eval(&z).unwrap();
                    rec.add_pt(*z, q_z).unwrap();
                }
            }
        }
    }}
    }
    panic!("Need more than 20 characteristics!");
}

fn rec_homogeneous_4(
    orig: &EvalCountingExpression<4>,
    mut rng: impl Rng,
) -> Rat<FlatPoly<Integer, 4>>
{
    use rare::rec::rat::cuyt_lee::{Rec4, Status};
    use rare::rec::rat::finite::cuyt_lee::Needed;
    const NVARS: usize = 4;

    let mut rec = Rec4::with_shift(1, rng.gen());
    seq!{ N in 0..20 {{
        const P: u64 = LARGE_PRIMES[N];
        let z: [Z64<P>; NVARS] = [(); NVARS].map(|_| rng.gen());
        let q_z = orig.try_eval(&z).unwrap();
        rec.add_pt(z, q_z).unwrap();
        loop {
            match rec.status().unwrap() {
                Status::Done => return rec.into_rat().unwrap(),
                Status::NextMod => break,
                Status::Needed(Needed::Pts(pts)) => {
                    let pts = Vec::from_iter(
                        pts.iter()
                            .filter_map(
                                |z| orig.try_eval(z).map(|q_z| (*z, q_z))
                            )
                    );
                    rec.add_pts(&pts).unwrap();
                },
                Status::Needed(Needed::Pt(z)) => {
                    let q_z: Z64<P> = orig.try_eval(&z).unwrap();
                    rec.add_pt(*z, q_z).unwrap();
                }
            }
        }
    }}
    }
    panic!("Need more than 20 characteristics!");
}

fn rec_scaling<const N: usize>(
    orig: &EvalCountingExpression<N>,
    mut rng: impl Rng,
) -> Rat<FlatPoly<Integer, N>> {
    use rare::rec::rat::thiele_multivar::{Rec, Status};
    use std::ops::ControlFlow::*;

    info!("Starting scaling Thiele reconstruction");
    const P: u64 = LARGE_PRIMES[0];

    let mut rec: Rec<P, N> = Rec::with_random_shift(2, &mut rng);
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
    seq!{ M in 0..20 {{
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
    panic!("Need more than 20 characteristics!");
}
