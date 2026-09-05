// Parse GenPept flat-file text into a compact TSV: version, locus_tags, old_locus_tags, genes, coded_by, definition
// Used both in the browser pane (NCBI efetch responses) and in Node (to verify against the Python parser).
function parseGenPept(t) {
  const out = [];
  for (const rec of t.split(/\n\/\/\n/)) {
    if (!/^LOCUS/m.test(rec)) continue;
    const ver = (rec.match(/^VERSION\s+(\S+)/m) || [])[1] || "";
    const lt = [...rec.matchAll(/\/locus_tag="([^"]+)"/g)].map(m => m[1]);
    const olt = [...rec.matchAll(/\/old_locus_tag="([^"]+)"/g)].map(m => m[1]);
    const genes = [...rec.matchAll(/\/gene="([^"]+)"/g)].map(m => m[1]);
    const cb = (rec.match(/\/coded_by="([^"]+)"/) || [])[1] || "";
    const defm = rec.match(/^DEFINITION\s+([\s\S]*?)\n(?=ACCESSION)/m);
    const def = defm ? defm[1].replace(/\n\s+/g, " ").trim() : "";
    out.push([ver, lt.join(";"), olt.join(";"), genes.join(";"), cb.replace(/\s+/g, ""), def].join("\t"));
  }
  return out.join("\n");
}
if (typeof module !== "undefined") {
  module.exports = { parseGenPept };
  if (require.main === module) {
    const fs = require("fs");
    for (const f of process.argv.slice(2)) process.stdout.write(parseGenPept(fs.readFileSync(f, "utf8")) + "\n");
  }
}
