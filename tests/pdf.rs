//! The hand-rolled PDF writer: structural validation (xref offsets,
//! stream lengths, page count) plus byte-exact determinism.

use kronrod::pdf::{render, render_dossier};

fn parse_xref_offsets(bytes: &[u8]) -> Vec<(usize, usize)> {
    // Returns (object_number, offset) for every "NNNNNNNNNN 00000 n" entry.
    let text = String::from_utf8_lossy(bytes);
    // Anchor on the standalone "xref\n" line (not "startxref\n").
    let xref_pos = text.rfind("\nxref\n").map(|p| p + 1).expect("xref table");
    let tail = &text[xref_pos..];
    let mut out = Vec::new();
    // rows after "xref\n" and the "0 {size}" subsection header: row i is
    // object i ("0 <offset> <gen> <f|n>").
    let mut obj = 0usize;
    for line in tail.lines().skip(2) {
        let line = line.trim_end();
        if line.is_empty() || line.starts_with("trailer") {
            break;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 3 && parts[2] == "n" {
            let offset: usize = parts[0].parse().unwrap();
            out.push((obj, offset));
        }
        obj += 1;
    }
    out
}

fn sample_report() -> String {
    let mut s = String::from("line one of the fake report\n");
    for i in 0..60 {
        s.push_str(&format!("body line {i:02} with some padding to be long\n"));
    }
    s
}

#[test]
fn pdf_header_and_footer() {
    let pdf = render("Test Title", "hello world");
    assert!(pdf.starts_with(b"%PDF-1.4\n"));
    let tail = &pdf[pdf.len() - 6..];
    assert_eq!(tail, b"%%EOF\n");
}

#[test]
fn xref_offsets_point_at_objects() {
    let pdf = render("Title (with parens) and backslash \\ text", &sample_report());
    let entries = parse_xref_offsets(&pdf);
    assert!(!entries.is_empty());
    for (obj, offset) in &entries {
        // "<obj> 0 obj" must start exactly at `offset`.
        let at = &pdf[*offset..];
        let expect = format!("{obj} 0 obj");
        assert!(
            at.starts_with(expect.as_bytes()),
            "object {obj} at {offset}: got {:?}",
            String::from_utf8_lossy(&at[..expect.len().min(at.len())])
        );
    }
}

#[test]
fn startxref_is_consistent() {
    let pdf = render("Title", "x");
    let text = String::from_utf8_lossy(&pdf).to_string();
    let pos = text.rfind("startxref\n").unwrap();
    let offset: usize = text[pos + 10..]
        .trim_end()
        .lines()
        .next()
        .unwrap()
        .parse()
        .unwrap();
    assert_eq!(&pdf[offset..offset + 5], b"xref\n");
}

#[test]
fn stream_lengths_match() {
    let pdf = render("Title", &sample_report());
    // Every "stream\n" block must be exactly /Length bytes before
    // "endstream".
    let bytes = &pdf[..];
    let mut pos = 0usize;
    let mut streams = 0usize;
    while let Some(rel) = bytes[pos..].windows(7).position(|w| w == b"stream\n") {
        let start = pos + rel + 7;
        let end_rel = bytes[start..]
            .windows(10)
            .position(|w| w == b"endstream\n")
            .expect("endstream");
        let len = end_rel;
        streams += 1;
        assert!(len > 0);
        // Extract /Length from the object header preceding the stream.
        let header = String::from_utf8_lossy(&bytes[..start]);
        let li = header.rfind("<< /Length ").unwrap();
        let len_decl: usize = header[li + 11..]
            .split(' ')
            .next()
            .unwrap()
            .parse()
            .unwrap();
        assert_eq!(len_decl, len, "stream {streams} length mismatch");
        pos = start + end_rel + 10; // skip past "endstream\n"
    }
    assert!(streams >= 1);
}

#[test]
fn page_count_and_kids() {
    let pdf = render("Title", &sample_report());
    let text = String::from_utf8_lossy(&pdf).to_string();
    // 60+ body lines > LINES_PER_PAGE → at least 2 pages.
    let count = text
        .split("/Count ")
        .nth(1)
        .and_then(|s| s.split('>').next())
        .and_then(|s| s.trim().parse::<usize>().ok())
        .expect("page count");
    let kids = text.matches("/Kids [").count();
    assert_eq!(kids, 1);
    assert!(count >= 2, "expected multiple pages, got {count}");
    // Each /Type /Page object references a /Contents stream.
    assert!(text.matches("/Type /Page ").count() == count);
}

#[test]
fn fixed_identity_no_clock() {
    let a = render_dossier(&sample_report());
    let b = render_dossier(&sample_report());
    assert_eq!(a, b, "dossier must be byte-identical");
    let text = String::from_utf8_lossy(&a).to_string();
    assert!(text.contains("D:20260101000000+00'00'"));
    assert!(text.contains("/ID [<4b524f4e442d444f53534945522d312e302e30>"));
}

#[test]
fn non_ascii_is_sanitised() {
    let pdf = render("Title π", "value = √2 ≈ 1.414 (π)");
    // No raw bytes ≥ 0x80 may appear in text content (the header comment
    // line is the only sanctioned binary region).
    let text = String::from_utf8_lossy(&pdf);
    // The binary comment is exactly one line right after the header.
    let after_header = &text[9..]; // skip "%PDF-1.4\n"
    let comment_end = after_header.find('\n').unwrap();
    let rest = &after_header[comment_end + 1..];
    assert!(rest.bytes().all(|b| b < 0x80), "non-ASCII leaked into PDF");
}

#[test]
fn title_with_parens_and_backslash_escapes() {
    let pdf = render("a (b) c \\ d", "x (y) z \\ w");
    let text = String::from_utf8_lossy(&pdf).to_string();
    assert!(text.contains("(a \\(b\\) c \\\\ d)"));
    // Balanced parentheses in the title string.
    let t = text
        .split("<< /Title (")
        .nth(1)
        .and_then(|s| s.split(") /Author").next())
        .unwrap();
    let open = t.chars().filter(|&c| c == '(').count();
    let close = t.chars().filter(|&c| c == ')').count();
    assert_eq!(open, close);
}
