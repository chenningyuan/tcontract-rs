use std::collections::BTreeSet;

use super::*;

#[test]
fn corpus_is_six_cases_at_three_aspect_ratios() {
    let cases = corpus();
    assert_eq!(cases.len(), 18);
    assert_eq!(
        cases.iter().map(|c| c.name).collect::<BTreeSet<_>>().len(),
        18
    );
    assert_eq!(
        cases.iter().map(|c| c.group).collect::<BTreeSet<_>>().len(),
        3
    );
    for group in ["ov_03", "ov_10", "ov_20"] {
        assert_eq!(cases.iter().filter(|c| c.group == group).count(), 6);
    }
    let expected = [
        "T2.vvvv",
        "T2.ring_mi",
        "T2.ring_mj",
        "T2.ring_u",
        "W.mbej.t2",
        "W.mbje.t2",
    ];
    for group in ["ov_03", "ov_10", "ov_20"] {
        for case in expected {
            assert!(
                cases.iter().any(|c| c.name == format!("{group}/{case}")),
                "missing {group}/{case}"
            );
        }
    }
}

#[test]
fn dfrccsd_einsums_match_the_paper_transcription() {
    let expected = [
        ("T2.vvvv", "ijef", "abef", "ijab"),
        ("T2.ring_mi", "miae", "mbje", "ijab"),
        ("T2.ring_mj", "mjae", "mbie", "ijab"),
        ("T2.ring_u", "imae", "mbej", "ijab"),
        ("W.mbej.t2", "jnbf", "mefn", "mbej"),
        ("W.mbje.t2", "njbf", "mefn", "mbje"),
    ];
    for (name, a, b, c) in expected {
        let full_name = format!("ov_03/{name}");
        let case = OV_03
            .iter()
            .find(|case| case.name == full_name)
            .unwrap_or_else(|| panic!("missing {full_name}"));
        assert_eq!((case.a, case.b, case.c), (a, b, c), "{full_name}");
    }
}

#[test]
fn ratio_extents_follow_the_derivation_rule() {
    // v = f·o with f the ratio denominator, and naux = round(k·(o+v)) at
    // the rung-median k of the NAO statistics panel.
    let expected = [
        ("ov_03", 3, 3.747, (18, 54, 270)),
        ("ov_10", 10, 2.568, (10, 100, 282)),
        ("ov_20", 20, 2.237, (7, 140, 329)),
    ];
    for (group, f, k, (o, v, naux)) in expected {
        let e = EXTENTS
            .iter()
            .find(|e| e.o == o && e.v == v && e.naux == naux)
            .unwrap_or_else(|| panic!("{group} extents missing from EXTENTS"));
        assert_eq!(e.v, f * e.o, "{group}: v must be f·o");
        assert_eq!(
            e.naux,
            (k * (e.o + e.v) as f64).round() as i64,
            "{group}: naux must follow the rung-median rule"
        );
        assert!(
            e.naux > e.v,
            "{group}: DF auxiliary space exceeds the virtual space"
        );
    }
}

#[test]
fn ov_03_extents_match_the_derived_numbers() {
    let cases = corpus();
    let t2 = size_case(&cases[0], 0.0);
    assert_eq!(t2.case.name, "ov_03/T2.vvvv");
    assert_eq!(t2.extent_of('i'), 18);
    assert_eq!(t2.extent_of('a'), 54);
    assert_eq!(t2.mnk(), (324, 2916, 2916));
    assert_eq!(t2.macs(), 2_754_990_144);
}

#[test]
fn layouts_and_mnk_agree_with_macs() {
    for case in corpus() {
        let sized = size_case(&case, 0.0);
        let (m, n, k) = sized.mnk();
        assert_eq!(m * n * k, sized.macs(), "case {}", case.name);
        let product = |indices: &str| {
            indices
                .chars()
                .map(|ch| sized.extent_of(ch) as usize)
                .product::<usize>()
        };
        assert_eq!(
            sized.elems_a(),
            product(case.a),
            "A layout of {}",
            case.name
        );
        assert_eq!(
            sized.elems_b(),
            product(case.b),
            "B layout of {}",
            case.name
        );
        assert_eq!(
            sized.elems_c(),
            product(case.c),
            "C layout of {}",
            case.name
        );
    }
}

#[test]
fn known_value_case_is_tiny_and_separate_from_the_corpus() {
    let known = known_value_case();
    assert_eq!(known.case.name, "known/T2.vvvv");
    assert_eq!(known.extent_of('a'), 3 * known.extent_of('i'));
    assert_eq!(known.mnk(), (4, 36, 36));
    assert_eq!(known.macs(), 5_184);
    assert_eq!(
        (known.elems_a(), known.elems_b(), known.elems_c()),
        (144, 1296, 144)
    );
    assert!(!corpus().iter().any(|case| case.name == known.case.name));
}

#[test]
fn whole_binary_footprints_stay_within_the_eight_gib_cell_budget() {
    // The extents are sized so the *whole binary* f64 footprint (A+B+C) of
    // every case fits an 8 GiB per-cell budget. The binding tensor of the
    // planned 25-case suite is the o^3 v^3 triples output of the
    // perturbative-triples cases, which this dfrccsd-only demo does not
    // carry -- so pin it explicitly: it is the tensor that consumes the
    // budget these extents were derived from, and it fits at every ratio.
    const CELL_BUDGET_ELEMS: i64 = 8 * 1024 * 1024 * 1024 / 8;
    for e in EXTENTS {
        let triples = e.o * e.o * e.o * e.v * e.v * e.v;
        assert!(
            triples <= CELL_BUDGET_ELEMS,
            "o^3 v^3 at (o={}, v={}) exceeds the 8 GiB cell budget",
            e.o,
            e.v
        );
    }
    for case in corpus() {
        let sized = size_case(&case, 0.0);
        let footprint = (sized.elems_a() + sized.elems_b() + sized.elems_c()) as i64;
        assert!(footprint <= CELL_BUDGET_ELEMS, "case {}", case.name);
    }
}
