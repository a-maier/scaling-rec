use std::collections::HashSet;

use anyhow::{Error, bail, Result};
use math_symbols::Symbol;

use crate::expression::{Atom, Expression, Operator};

impl<'a, const N: usize> TryFrom<&'a str> for Expression<N> {
    type Error = anyhow::Error;

    fn try_from(s: &'a str) -> Result<Self, Self::Error> {
        let mut parser = Parser::default();
        parser.parse(s)?;
        parser.try_into()
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct Parser {
    output: Vec<Atom>,
    operators: Vec<Operator>,
    symbols: HashSet<Symbol>,
}

impl Parser {
    //shunting-yard parser
    fn parse(&mut self, mut s: &str) -> Result<()> {
        s = s.trim_start();
        while !s.is_empty() {
            if s.starts_with(|c: char| c.is_ascii_digit()) {
                let int_end = s.find(|c: char| !c.is_ascii_digit())
                    .unwrap_or(s.len());
                let int;
                (int, s) = s.split_at(int_end);
                let int = int.parse()?;
                self.output.push(Atom::Number(int));
            } else if s.starts_with(|c: char| c.is_alphabetic()) {
                let sym_end = s.find(|c: char| c != '_' && !c.is_alphanumeric())
                    .unwrap_or(s.len());
                let sym;
                (sym, s) = s.split_at(sym_end);
                let sym = Symbol::new(sym);
                self.symbols.insert(sym);
                self.output.push(sym.into());
            } else {
                use Operator::*;
                match s.bytes().next().unwrap() {
                    b'+' => self.push_op(Plus),
                    b'-' => self.push_op(Minus),
                    b'*' => self.push_op(Times),
                    b'/' => self.push_op(Divided),
                    b'^' => self.push_op(Power),
                    b'(' => self.operators.push(LeftBracket),
                    b')' => self.close_bracket()?,
                    _ => bail!("Failed to parse start of {s}")
                }
                s = &s[1..]
            }
            s = s.trim_start();
        }
        Ok(())
    }

    fn push_op(&mut self, op: Operator) {
        let prec = op.precedence();
        // all operators are treated as left-associative
        // and left parentheses have the lowest precedence
        // so this check is enough
        while self.last_op_precedence() >=  prec {
            let op = self.operators.pop().unwrap();
            self.output.push(op.into());
        }
        self.operators.push(op);
    }

    fn last_op_precedence(&self) -> u32 {
        self.operators.last()
            .map(|op| op.precedence())
            .unwrap_or(0)
    }

    fn close_bracket(&mut self) -> Result<()> {
        loop {
            let Some(op) = self.operators.pop() else {
                bail!("Closing bracket without matching opening bracket")
            };
            if op == Operator::LeftBracket {
                return Ok(());
            }
            self.output.push(op.into())
        }
    }
}

impl<const N: usize> TryFrom<Parser> for Expression<N> {
    type Error = Error;

    fn try_from(parser: Parser) -> std::result::Result<Self, Self::Error> {
        let Parser { mut output, mut operators, symbols } = parser;
        let mut symbols = Vec::from_iter(symbols);
        symbols.sort();
        while let Some(op) = operators.pop() {
            if op == Operator::LeftBracket {
                bail!("Unclosed opening bracket")
            }
            output.push(op.into())
        }
        let symbols = symbols.try_into().unwrap();
        Ok(Expression{
            atoms: output,
            symbols
        })
    }
}
