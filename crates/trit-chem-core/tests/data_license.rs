//! Licence and provenance gates for the element table.
//!
//! The first element source was CC BY-SA 3.0 (ShareAlike would propagate into
//! this MIT repository). These tests make the replacement permanent: the MIT
//! notice must be present, the derived JSON must declare MIT, and the old
//! source file must not come back.

use std::path::Path;

fn root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

#[test]
fn mendeleev_mit_notice_is_present() {
    let licence = std::fs::read_to_string(root().join("data/external/MENDELEEV-LICENSE"))
        .expect("MENDELEEV-LICENSE must be committed");
    assert!(licence.contains("MIT License"), "not an MIT notice");
    assert!(
        licence.contains("Lukasz Mentel"),
        "copyright holder missing from the notice"
    );
}

#[test]
fn derived_element_json_declares_mit_and_pins_the_wheel() {
    let json = std::fs::read_to_string(root().join("data/external/elements.json"))
        .expect("elements.json must be committed");
    assert!(json.contains("\"license\": \"MIT\""));
    assert!(json.contains("\"package\": \"mendeleev\""));
    // 64 hex characters: the SHA-256 of the pinned wheel.
    let hash = json
        .split("\"wheel_sha256\": \"")
        .nth(1)
        .and_then(|s| s.split('"').next())
        .expect("wheel hash");
    assert_eq!(hash.len(), 64, "hash is not SHA-256");
}

#[test]
fn share_alike_source_is_gone() {
    assert!(
        !root().join("data/external/PeriodicTableJSON.json").exists(),
        "CC BY-SA element source is back; do not vendor it"
    );
}

#[test]
fn generated_table_names_its_source() {
    let src = std::fs::read_to_string(root().join("crates/trit-chem-core/src/elements_data.rs"))
        .expect("generated table");
    assert!(
        src.contains("mendeleev"),
        "generated header lost its source"
    );
    assert!(!src.contains("PeriodicTableJSON"));
}
