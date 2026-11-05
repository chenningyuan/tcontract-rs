//! Six leading DF-RCCSD binary contractions at fixed QC aspect ratios.
//!
//! Equations and their order come from U. Bozkaya and C. D. Sherrill,
//! J. Chem. Phys. 147, 044104 (2017), DOI 10.1063/1.4994918, main Eq. (5)
//! and supplementary Eqs. (10)-(43). The six cases are the two leading
//! formal-scaling rungs of the DF-RCCSD energy/amplitude route: the
//! `o^2 v^4` head (`T2.vvvv`) and the five `o^3 v^3` ring terms. Occupied
//! indices use `o`, virtual indices `v`, and the density-fitting index `Q`
//! uses `naux` (no case in this demo carries `Q`; the extents keep the DF
//! domain so the same triples serve the rest of the planned suite).
//!
//! Extents are aspect-ratio classes, not molecules. The three
//! occupied/virtual ratios o/v = 1/3, 1/10, 1/20 are the median v/o of the
//! cc-pVDZ, cc-pVTZ, and cc-pVQZ rungs of a 5,239-molecule NAO statistics
//! panel (3.90, 10.17, 20.36), so each ratio stands for one basis-quality
//! rung without pinning a molecule, and `naux = round(k_rung*(o+v))` with
//! k_rung the rung median of naux/nao (3.747, 2.568, 2.237).
//!
//! The extents are derived from a whole-binary memory budget: the largest
//! f64 A+B+C footprint over the planned 25-case suite this demo belongs to
//! -- bound by the `o^3 v^3` triples output of the perturbative-triples
//! cases -- stays within 8 GiB per cell, so no case is skipped for memory.
//! At that budget: (o, v, naux) = (18, 54, 270), (10, 100, 282),
//! (7, 140, 329), and `v*o` stays nearly constant (972 / 1000 / 980): a
//! larger basis rung affords more virtuals but fewer occupied orbitals.

#[derive(Clone, Copy, Debug)]
pub struct Extents {
    pub o: i64,
    pub v: i64,
    pub naux: i64,
}

static EXTENTS: &[Extents] = &[
    // o/v = 1/3, the cc-pVDZ rung median (v/o 3.90, naux/nao 3.747).
    Extents {
        o: 18,
        v: 54,
        naux: 270,
    },
    // o/v = 1/10, the cc-pVTZ rung median (v/o 10.17, naux/nao 2.568).
    Extents {
        o: 10,
        v: 100,
        naux: 282,
    },
    // o/v = 1/20, the cc-pVQZ rung median (v/o 20.36, naux/nao 2.237).
    Extents {
        o: 7,
        v: 140,
        naux: 329,
    },
];

#[derive(Clone, Debug)]
pub struct Layout {
    extents: Vec<i64>,
    strides: Vec<i64>,
}

impl Layout {
    pub fn col_major(extents: &[i64]) -> Self {
        let mut strides = Vec::with_capacity(extents.len());
        let mut acc = 1i64;
        for &e in extents {
            strides.push(acc);
            acc *= e.max(1);
        }
        Self {
            extents: extents.to_vec(),
            strides,
        }
    }

    pub fn extents(&self) -> &[i64] {
        &self.extents
    }

    pub fn strides(&self) -> &[i64] {
        &self.strides
    }

    pub fn storage_len(&self) -> i64 {
        let mut n = 1i64;
        for (&e, &s) in self.extents.iter().zip(&self.strides) {
            if e > 0 {
                n += (e - 1) * s.abs();
            }
        }
        n
    }
}

#[derive(Clone, Debug)]
pub struct Case {
    pub name: &'static str,
    pub group: &'static str,
    pub c: &'static str,
    pub a: &'static str,
    pub b: &'static str,
    pub extents: &'static Extents,
}

// Dedicated tiny extents for the shared verifier's all-ones check, keeping
// the ov_03 aspect ratio o/v = 1/3 and its naux rule:
// round(3.747 * (2 + 6)) = 30. This case is outside the corpus and never
// changes the fixed QC extents above.
static KNOWN_EXTENTS: Extents = Extents {
    o: 2,
    v: 6,
    naux: 30,
};
static KNOWN_CASE: Case = Case {
    name: "known/T2.vvvv",
    group: "known",
    c: "ijab",
    a: "ijef",
    b: "abef",
    extents: &KNOWN_EXTENTS,
};

#[derive(Clone, Debug)]
pub struct Sized {
    pub case: Case,
    pub labels: Vec<char>,
    pub extents: Vec<i64>,
    pub la: Layout,
    pub lb: Layout,
    pub lc: Layout,
    pub idx_a: Vec<i64>,
    pub idx_b: Vec<i64>,
    pub idx_c: Vec<i64>,
}

impl Sized {
    #[cfg(test)]
    pub fn extent_of(&self, ch: char) -> i64 {
        self.extents[self.labels.iter().position(|&c| c == ch).unwrap()]
    }

    pub fn macs(&self) -> u64 {
        self.extents.iter().map(|&e| e as u64).product()
    }

    pub fn elems_a(&self) -> usize {
        self.la.storage_len() as usize
    }

    pub fn elems_b(&self) -> usize {
        self.lb.storage_len() as usize
    }

    pub fn elems_c(&self) -> usize {
        self.lc.storage_len() as usize
    }

    pub fn mnk(&self) -> (u64, u64, u64) {
        let inm = |ch: char| self.case.c.contains(ch) && self.case.a.contains(ch);
        let inn = |ch: char| self.case.c.contains(ch) && self.case.b.contains(ch);
        let mut m = 1u64;
        let mut n = 1u64;
        let mut k = 1u64;
        for (i, &ch) in self.labels.iter().enumerate() {
            let e = self.extents[i] as u64;
            if inm(ch) {
                m *= e;
            } else if inn(ch) {
                n *= e;
            } else {
                k *= e;
            }
        }
        (m, n, k)
    }
}

macro_rules! qc_cases {
    ($extents:expr, $group:literal) => {
        &[
            Case {
                name: concat!($group, "/T2.vvvv"),
                group: $group,
                c: "ijab",
                a: "ijef",
                b: "abef",
                extents: &EXTENTS[$extents],
            },
            Case {
                name: concat!($group, "/T2.ring_mi"),
                group: $group,
                c: "ijab",
                a: "miae",
                b: "mbje",
                extents: &EXTENTS[$extents],
            },
            Case {
                name: concat!($group, "/T2.ring_mj"),
                group: $group,
                c: "ijab",
                a: "mjae",
                b: "mbie",
                extents: &EXTENTS[$extents],
            },
            Case {
                name: concat!($group, "/T2.ring_u"),
                group: $group,
                c: "ijab",
                a: "imae",
                b: "mbej",
                extents: &EXTENTS[$extents],
            },
            Case {
                name: concat!($group, "/W.mbej.t2"),
                group: $group,
                c: "mbej",
                a: "jnbf",
                b: "mefn",
                extents: &EXTENTS[$extents],
            },
            Case {
                name: concat!($group, "/W.mbje.t2"),
                group: $group,
                c: "mbje",
                a: "njbf",
                b: "mefn",
                extents: &EXTENTS[$extents],
            },
        ]
    };
}

static OV_03: &[Case] = qc_cases!(0, "ov_03");
static OV_10: &[Case] = qc_cases!(1, "ov_10");
static OV_20: &[Case] = qc_cases!(2, "ov_20");

pub fn corpus() -> Vec<Case> {
    [OV_03, OV_10, OV_20]
        .into_iter()
        .flat_map(|cases| cases.iter().cloned())
        .collect()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Stress {
    #[default]
    None,
}

impl Stress {
    pub fn parse(s: &str) -> Option<Stress> {
        (s == "none").then_some(Stress::None)
    }

    pub fn name(self) -> &'static str {
        "none"
    }
}

/// The tiny all-ones contraction used by the shared verifier. Normal corpus
/// cases always use their declared extents and never consult this helper.
pub fn known_value_case() -> Sized {
    size_case(&KNOWN_CASE, 0.0)
}

pub fn size_case(case: &Case, tensor_bytes: f64) -> Sized {
    size_case_stressed(case, tensor_bytes, Stress::None)
}

pub fn size_case_stressed(case: &Case, tensor_bytes: f64, stress: Stress) -> Sized {
    let _ = tensor_bytes;
    let Stress::None = stress;
    let extent = |ch: char| match ch {
        'i' | 'j' | 'k' | 'l' | 'm' | 'n' => case.extents.o,
        'a' | 'b' | 'c' | 'd' | 'e' | 'f' => case.extents.v,
        'Q' => case.extents.naux,
        _ => panic!("unknown QC index {ch} in {}", case.name),
    };
    let mut labels = Vec::new();
    for indices in [case.c, case.a, case.b] {
        for ch in indices.chars() {
            if !labels.contains(&ch) {
                labels.push(ch);
            }
        }
    }
    let extents: Vec<i64> = labels.iter().map(|&ch| extent(ch)).collect();
    let shape = |indices: &str| -> Vec<i64> { indices.chars().map(extent).collect() };
    let idx = |indices: &str| -> Vec<i64> { indices.chars().map(|ch| ch as i64).collect() };
    Sized {
        case: case.clone(),
        labels,
        extents,
        la: Layout::col_major(&shape(case.a)),
        lb: Layout::col_major(&shape(case.b)),
        lc: Layout::col_major(&shape(case.c)),
        idx_a: idx(case.a),
        idx_b: idx(case.b),
        idx_c: idx(case.c),
    }
}

#[cfg(test)]
mod tests;
