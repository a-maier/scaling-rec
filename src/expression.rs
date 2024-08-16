use std::cell::RefCell;

use log::info;
use math_symbols::Symbol;
use rare::{Z64, traits::{TryEval, Zero}};
use rug::{Integer, integer::IntegerExt64};

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Atom {
    Symbol(Symbol),
    Number(Integer),
    Operator(Operator),
}

impl From<Operator> for Atom {
    fn from(op: Operator) -> Self {
        Self::Operator(op)
    }
}

impl From<Integer> for Atom {
    fn from(int: Integer) -> Self {
        Self::Number(int)
    }
}

impl From<Symbol> for Atom {
    fn from(sym: Symbol) -> Self {
        Self::Symbol(sym)
    }
}


#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Operator {
    Plus,
    Minus,
    Times,
    Divided,
    Power,
    LeftBracket,
    // unused: RightBracket,
}

impl Operator {
    pub(crate) fn precedence(&self) -> u32 {
        use Operator::*;
        match self {
            Plus => 10,
            Minus => 10,
            Times => 20,
            Divided => 20,
            Power => 30,
            LeftBracket => 0,
        }
    }
}

// expression counting the number of evaluations
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct EvalCountingExpression<const N: usize>{
    last_mod: RefCell<u64>,
    pub(crate) count: RefCell<usize>,
    expression: Expression<N>,
}

impl<const N: usize> EvalCountingExpression<N> {
    pub(crate) fn print_and_reset_count(&self) {
        let Self { last_mod, count, expression: _ } = self;
        if *count.borrow() != 0 {
            info!("Probes mod {}: {}", last_mod.borrow(), count.borrow())
        }
        *count.borrow_mut() = 0;
    }
}

impl<'a, const N: usize> TryFrom<&'a str> for EvalCountingExpression<N> {
    type Error = anyhow::Error;

    fn try_from(s: &'a str) -> Result<Self, Self::Error> {
        let expression = Expression::try_from(s)?;
        Ok(Self {
            last_mod: RefCell::new(0),
            expression,
            count: RefCell::new(0),
        })
    }
}

impl<const P: u64, const N: usize> TryEval<[Z64<P>; N]> for EvalCountingExpression<N> {
    type Output = Z64<P>;

    fn try_eval(&self, pt: &[Z64<P>; N]) -> Option<Self::Output> {
        let Self { last_mod, count, expression } = self;
        if P != *last_mod.borrow() {
            self.print_and_reset_count();
            *last_mod.borrow_mut() = P;
        }
        *count.borrow_mut() += 1;
        expression.try_eval(pt)
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Expression<const N: usize> {
    pub(crate) symbols: [Symbol; N],
    pub(crate) atoms: Vec<Atom>,
}

impl<const P: u64, const N: usize> TryEval<[Z64<P>; N]> for Expression<N> {
    type Output = Z64<P>;

    fn try_eval(&self, pt: &[Z64<P>; N]) -> Option<Self::Output> {
        let mut stack = Vec::new();
        for atom in &self.atoms {
            match atom {
                Atom::Symbol(s) => {
                    let idx = self.symbols.iter().position(|x| x == s).unwrap();
                    stack.push(pt[idx])
                },
                Atom::Number(n) => {
                    let n = unsafe { Z64::new_unchecked(n.mod_u64(P)) };
                    stack.push(n)
                },
                Atom::Operator(op) => {
                    let op = *op;
                    use Operator::*;
                    let Some(op2) = stack.pop() else {
                        panic!("No operand for {op:?} on stack!");
                    };
                    let op1 = match stack.pop() {
                        Some(op) => op,
                        None => if op == Minus || op == Plus {
                            Z64::zero()
                        } else {
                            panic!("No operand for {op:?} on stack!");
                        },
                    };
                    let res = match op {
                        Plus => op1 + op2,
                        Minus => op1 - op2,
                        Times => op1 * op2,
                        Divided => {
                            let inv_op2 = op2.try_inv()?;
                            op1 * inv_op2
                        },
                        Power => op1.powu(op2.into()),
                        LeftBracket => unreachable!("Opening bracket in expression"),
                    };
                    stack.push(res);
                },
            }
        }
        if stack.len() == 1 {
            Some(stack[0])
        } else {
            panic!("Stack has {} elements left after eval", stack.len());
        }
    }
}
