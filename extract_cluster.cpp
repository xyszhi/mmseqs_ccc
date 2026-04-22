#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <sys/stat.h>

// ---------------------------------------------------------------------------
// Usage
// ---------------------------------------------------------------------------
static void print_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " [OPTIONS] <rep_name> [rep_name2 ...]\n"
        "\n"
        "Extract alignment results for all members of specified cluster(s).\n"
        "\n"
        "Required arguments:\n"
        "  <rep_name>          Representative sequence name(s) (column 1 of clusters TSV)\n"
        "\n"
        "Options:\n"
        "  --db        <path>  DB prefix (default: test/filtered_db)\n"
        "  --lookup    <path>  Lookup file (default: test/protein_db.lookup)\n"
        "  --clusters  <path>  Clusters TSV (default: <db>_clusters.tsv)\n"
        "  --format    <fmt>   Output format string (default: query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,qcov,tstart,tend,tlen,tcov,evalue,bits)\n"
        "                      Available fields: query target fident alnlen evalue\n"
        "                                        qstart qend qlen tstart tend tlen\n"
        "                                        qcov tcov mismatch gapopen bits\n"
        "                      Note: alnlen/mismatch/gapopen require DB built with 'mmseqs search -a'.\n"
        "                            mismatch is approximate (fident stored at 3 decimal precision).\n"
        "  --sep       <char>  Output field separator (default: tab)\n"
        "  --out       <path>  Output file (default: stdout)\n"
        "  --header            Print header line\n"
        "  --help              Show this help\n"
        "\n"
        "Format string example: --format \"query,target,fident,evalue,qcov,tcov\"\n";
}

// ---------------------------------------------------------------------------
// Fast integer parser
// ---------------------------------------------------------------------------
static inline const char* parse_uint(const char* p, const char* end, long long& out) {
    while (p < end && (*p == ' ' || *p == '\t')) ++p;
    if (p >= end || *p < '0' || *p > '9') return nullptr;
    long long v = 0;
    while (p < end && *p >= '0' && *p <= '9') { v = v * 10 + (*p - '0'); ++p; }
    out = v;
    return p;
}

// ---------------------------------------------------------------------------
// Index entry
// ---------------------------------------------------------------------------
struct IndexEntry {
    int       query_id;
    long long local_offset;
    long long length;
};

// ---------------------------------------------------------------------------
// Parse cigar string to compute alnlen and gapopen
// cigar format: e.g. "88M3D242M1I10M"
// M = match/mismatch column, D = deletion (gap in query), I = insertion (gap in target)
// alnlen  = sum of all M+D+I lengths
// gapopen = number of D runs + number of I runs
// ---------------------------------------------------------------------------
static void parse_cigar(const char* s, const char* e, int& alnlen, int& gapopen) {
    alnlen  = 0;
    gapopen = 0;
    int num = 0;
    for (const char* p = s; p < e; ++p) {
        char c = *p;
        if (c >= '0' && c <= '9') {
            num = num * 10 + (c - '0');
        } else {
            if (num == 0) num = 1;
            if (c == 'M' || c == 'D' || c == 'I') {
                alnlen += num;
                if (c == 'D' || c == 'I') gapopen++;
            }
            num = 0;
        }
    }
}

// ---------------------------------------------------------------------------
// Alignment record (raw fields from DB)
// Internal mmseqs2 format: target_id score fident evalue qstart qend qlen tstart tend tlen [cigar]
// score = bit score (as computed by evaluer->computeBitScore)
// ---------------------------------------------------------------------------
struct AlnRecord {
    int    query_id;
    int    target_id;
    int    raw_score;  // bit score (mmseqs2 evaluer->computeBitScore result)
    double fident;
    double evalue;
    int    qstart;
    int    qend;
    int    qlen;
    int    tstart;
    int    tend;
    int    tlen;
    // derived from cigar (if available, else -1)
    int    alnlen;
    int    gapopen;
    int    mismatch;
};

// ---------------------------------------------------------------------------
// Format a single record according to format fields
// ---------------------------------------------------------------------------
static void write_record(std::FILE* out,
                         const AlnRecord& r,
                         const std::string& qname,
                         const std::string& tname,
                         const std::vector<std::string>& fields,
                         const std::string& sep)
{
    for (size_t i = 0; i < fields.size(); ++i) {
        if (i > 0) std::fputs(sep.c_str(), out);
        const std::string& f = fields[i];
        if (f == "query")    std::fputs(qname.c_str(), out);
        else if (f == "target")   std::fputs(tname.c_str(), out);
        else if (f == "fident")   std::fprintf(out, "%.3f", r.fident);
        else if (f == "alnlen")   { if (r.alnlen >= 0) std::fprintf(out, "%d", r.alnlen); else std::fputs("N/A", out); }
        else if (f == "evalue")   std::fprintf(out, "%g", r.evalue);
        else if (f == "qstart")   std::fprintf(out, "%d", r.qstart + 1);
        else if (f == "qend")     std::fprintf(out, "%d", r.qend + 1);
        else if (f == "qlen")     std::fprintf(out, "%d", r.qlen);
        else if (f == "tstart")   std::fprintf(out, "%d", r.tstart + 1);
        else if (f == "tend")     std::fprintf(out, "%d", r.tend + 1);
        else if (f == "tlen")     std::fprintf(out, "%d", r.tlen);
        else if (f == "qcov") {
            double qcov = r.qlen > 0 ? (double)(r.qend - r.qstart) / r.qlen : 0.0;
            std::fprintf(out, "%.4f", qcov);
        }
        else if (f == "tcov") {
            double tcov = r.tlen > 0 ? (double)(r.tend - r.tstart) / r.tlen : 0.0;
            std::fprintf(out, "%.4f", tcov);
        }
        else if (f == "mismatch") { if (r.mismatch >= 0) std::fprintf(out, "%d", r.mismatch); else std::fputs("N/A", out); }
        else if (f == "gapopen")  { if (r.gapopen >= 0) std::fprintf(out, "%d", r.gapopen); else std::fputs("N/A", out); }
        else if (f == "bits")     std::fprintf(out, "%d", r.raw_score);
        else std::fputs(f.c_str(), out); // unknown field: output as-is
    }
    std::fputc('\n', out);
}

// ---------------------------------------------------------------------------
// Split format string by comma
// ---------------------------------------------------------------------------
static std::vector<std::string> split_format(const std::string& fmt) {
    std::vector<std::string> result;
    std::stringstream ss(fmt);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        // trim whitespace
        size_t a = tok.find_first_not_of(" \t");
        size_t b = tok.find_last_not_of(" \t");
        if (a != std::string::npos)
            result.push_back(tok.substr(a, b - a + 1));
    }
    return result;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    std::string db_prefix    = "test/filtered_db";
    std::string lookup_path  = "";
    std::string clusters_path = "";
    std::string format_str   = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,qcov,tstart,tend,tlen,tcov,evalue,bits";
    std::string sep          = "\t";
    std::string out_path     = "";
    bool print_header        = false;
    std::vector<std::string> rep_names;

    // Parse arguments
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(argv[0]); return 0; }
        else if (a == "--db"       && i+1 < argc) { db_prefix    = argv[++i]; }
        else if (a == "--lookup"   && i+1 < argc) { lookup_path  = argv[++i]; }
        else if (a == "--clusters" && i+1 < argc) { clusters_path = argv[++i]; }
        else if (a == "--format"   && i+1 < argc) { format_str   = argv[++i]; }
        else if (a == "--sep"      && i+1 < argc) { sep          = argv[++i]; }
        else if (a == "--out"      && i+1 < argc) { out_path     = argv[++i]; }
        else if (a == "--header")                  { print_header = true; }
        else if (a[0] != '-')                      { rep_names.push_back(a); }
        else { std::cerr << "Unknown option: " << a << "\n"; print_usage(argv[0]); return 1; }
    }

    if (rep_names.empty()) {
        std::cerr << "Error: no representative sequence name(s) specified.\n\n";
        print_usage(argv[0]);
        return 1;
    }

    if (lookup_path.empty())   lookup_path   = db_prefix.substr(0, db_prefix.rfind('/') + 1) + "protein_db.lookup";
    if (clusters_path.empty()) clusters_path = db_prefix + "_clusters.tsv";

    // Handle \t in sep
    if (sep == "\\t") sep = "\t";

    std::vector<std::string> fields = split_format(format_str);
    if (fields.empty()) { std::cerr << "Error: empty format string.\n"; return 1; }

    // -----------------------------------------------------------------------
    // Load lookup: name -> id  and  id -> name
    // -----------------------------------------------------------------------
    std::unordered_map<std::string, int> name2id;
    std::unordered_map<int, std::string> id2name;
    {
        std::ifstream lf(lookup_path);
        if (!lf.is_open()) {
            std::cerr << "Error: cannot open lookup file: " << lookup_path << "\n";
            return 1;
        }
        std::string line;
        while (std::getline(lf, line)) {
            if (line.empty()) continue;
            const char* p = line.data();
            const char* end = p + line.size();
            long long id;
            p = parse_uint(p, end, id);
            if (!p) continue;
            while (p < end && (*p == ' ' || *p == '\t')) ++p;
            const char* ns = p;
            while (p < end && *p != '\r' && *p != '\n' && *p != '\t') ++p;
            std::string name(ns, p);
            name2id[name] = (int)id;
            id2name[(int)id] = name;
        }
    }
    std::cerr << "Loaded " << name2id.size() << " sequences from lookup.\n";

    // -----------------------------------------------------------------------
    // Load clusters TSV: find all members for requested rep names
    // -----------------------------------------------------------------------
    std::unordered_set<std::string> rep_set(rep_names.begin(), rep_names.end());
    // query_id -> set of target_ids to extract
    // We collect all (query_id, target_id) pairs we want
    // Actually: for each member sequence, we need to find its alignment block (query_id = member_id)
    // The DB stores: for each query, all its alignment targets
    // So we need to collect all member IDs, then for each member find its block in the DB

    // member_ids: set of query IDs whose blocks we want to extract
    std::unordered_set<int> member_ids;
    // Also track which rep each member belongs to (for display if needed)
    std::unordered_map<int, std::string> member_to_rep; // member_id -> rep_name

    {
        std::ifstream cf(clusters_path);
        if (!cf.is_open()) {
            std::cerr << "Error: cannot open clusters file: " << clusters_path << "\n";
            return 1;
        }
        std::string line;
        while (std::getline(cf, line)) {
            if (line.empty()) continue;
            size_t tab = line.find('\t');
            if (tab == std::string::npos) continue;
            std::string rep    = line.substr(0, tab);
            std::string member = line.substr(tab + 1);
            // trim \r
            if (!member.empty() && member.back() == '\r') member.pop_back();
            if (rep_set.count(rep)) {
                auto it = name2id.find(member);
                if (it != name2id.end()) {
                    member_ids.insert(it->second);
                    member_to_rep[it->second] = rep;
                } else {
                    std::cerr << "Warning: member '" << member << "' not found in lookup.\n";
                }
            }
        }
    }

    if (member_ids.empty()) {
        std::cerr << "Error: no members found for the specified representative(s).\n";
        return 1;
    }
    std::cerr << "Found " << member_ids.size() << " member sequences to extract.\n";

    // -----------------------------------------------------------------------
    // Read index
    // -----------------------------------------------------------------------
    int num_splits = 0;
    {
        struct stat st;
        while (stat((db_prefix + "." + std::to_string(num_splits)).c_str(), &st) == 0)
            num_splits++;
    }
    if (num_splits == 0) {
        std::cerr << "Error: no split files found for prefix: " << db_prefix << "\n";
        return 1;
    }

    std::vector<long long> cum_sizes(num_splits + 1, 0);
    for (int i = 0; i < num_splits; i++) {
        struct stat st;
        stat((db_prefix + "." + std::to_string(i)).c_str(), &st);
        cum_sizes[i + 1] = cum_sizes[i] + (long long)st.st_size;
    }

    // file_id -> list of IndexEntry for members we want
    std::vector<std::vector<IndexEntry>> file_entries(num_splits);
    {
        std::ifstream idx(db_prefix + ".index");
        if (!idx.is_open()) {
            std::cerr << "Error: cannot open index: " << db_prefix << ".index\n";
            return 1;
        }
        std::string line;
        while (std::getline(idx, line)) {
            if (line.empty()) continue;
            const char* p = line.data();
            const char* end = p + line.size();
            long long qid, offset, length;
            p = parse_uint(p, end, qid);
            if (!p) continue;
            while (p < end && (*p == ' ' || *p == '\t')) ++p;
            const char* p2 = parse_uint(p, end, offset);
            if (!p2) continue; p = p2;
            while (p < end && (*p == ' ' || *p == '\t')) ++p;
            p2 = parse_uint(p, end, length);
            if (!p2) continue;

            if (!member_ids.count((int)qid)) continue;

            int fid = (int)(std::upper_bound(cum_sizes.begin(), cum_sizes.end(), offset)
                            - cum_sizes.begin()) - 1;
            if (fid < 0) fid = 0;
            if (fid >= num_splits) fid = num_splits - 1;
            file_entries[fid].push_back({(int)qid, offset - cum_sizes[fid], length});
        }
    }

    // -----------------------------------------------------------------------
    // Open output
    // -----------------------------------------------------------------------
    std::FILE* fout = stdout;
    if (!out_path.empty()) {
        fout = std::fopen(out_path.c_str(), "w");
        if (!fout) { std::cerr << "Error: cannot open output file: " << out_path << "\n"; return 1; }
    }

    // Print header
    if (print_header) {
        for (size_t i = 0; i < fields.size(); ++i) {
            if (i > 0) std::fputs(sep.c_str(), fout);
            std::fputs(fields[i].c_str(), fout);
        }
        std::fputc('\n', fout);
    }

    // -----------------------------------------------------------------------
    // Extract and output records
    // -----------------------------------------------------------------------
    long long total_records = 0;
    const int IO_BUF = 16 << 20;
    std::vector<char> io_buf(IO_BUF);

    for (int fid = 0; fid < num_splits; ++fid) {
        auto& entries = file_entries[fid];
        if (entries.empty()) continue;

        std::sort(entries.begin(), entries.end(),
                  [](const IndexEntry& a, const IndexEntry& b) {
                      return a.local_offset < b.local_offset;
                  });

        std::ifstream sf(db_prefix + "." + std::to_string(fid), std::ios::binary);
        if (!sf.is_open()) {
            std::cerr << "Warning: cannot open split file " << fid << "\n";
            continue;
        }

        long long buf_start = -1, buf_end = -1;

        auto ensure_buf = [&](long long need_start, long long need_len) -> bool {
            if (need_start >= buf_start && need_start + need_len <= buf_end) return true;
            sf.seekg(need_start);
            if (!sf) return false;
            long long to_read = std::max(need_len, (long long)IO_BUF);
            sf.read(io_buf.data(), to_read);
            long long got = sf.gcount();
            if (got <= 0) return false;
            buf_start = need_start;
            buf_end   = need_start + got;
            return (got >= need_len);
        };

        for (auto& e : entries) {
            const char* block_ptr = nullptr;
            std::string fallback_buf;

            if (e.length <= IO_BUF) {
                if (!ensure_buf(e.local_offset, e.length)) continue;
                block_ptr = io_buf.data() + (e.local_offset - buf_start);
            } else {
                fallback_buf.resize(e.length);
                sf.seekg(e.local_offset);
                sf.read(&fallback_buf[0], e.length);
                if (sf.gcount() < e.length) continue;
                block_ptr = fallback_buf.data();
                buf_start = buf_end = -1;
            }

            const std::string& qname = id2name.count(e.query_id) ? id2name[e.query_id] : std::to_string(e.query_id);

            const char* p   = block_ptr;
            const char* end = block_ptr + e.length;

            while (p < end) {
                if (*p == '\0') { ++p; continue; }
                const char* lend = p;
                while (lend < end && *lend != '\n') ++lend;

                // Parse 10 or 11 tab-separated fields:
                // target_id raw_score fident evalue qstart qend qlen tstart tend tlen [cigar]
                auto next_tab = [](const char* s, const char* e) -> const char* {
                    while (s < e && *s != '\t') ++s;
                    return s;
                };

                const char* f1 = p;
                const char* f1e = next_tab(f1, lend);
                const char* f2 = f1e < lend ? f1e + 1 : lend;
                const char* f2e = next_tab(f2, lend);
                const char* f3 = f2e < lend ? f2e + 1 : lend;
                const char* f3e = next_tab(f3, lend);
                const char* f4 = f3e < lend ? f3e + 1 : lend;
                const char* f4e = next_tab(f4, lend);
                const char* f5 = f4e < lend ? f4e + 1 : lend;
                const char* f5e = next_tab(f5, lend);
                const char* f6 = f5e < lend ? f5e + 1 : lend;
                const char* f6e = next_tab(f6, lend);
                const char* f7 = f6e < lend ? f6e + 1 : lend;
                const char* f7e = next_tab(f7, lend);
                const char* f8 = f7e < lend ? f7e + 1 : lend;
                const char* f8e = next_tab(f8, lend);
                const char* f9 = f8e < lend ? f8e + 1 : lend;
                const char* f9e = next_tab(f9, lend);
                const char* f10 = f9e < lend ? f9e + 1 : lend;
                const char* f10e = next_tab(f10, lend);
                // optional 11th field: cigar
                const char* f11 = f10e < lend ? f10e + 1 : lend;
                const char* f11e = next_tab(f11, lend);

                if (f1e > f1) {
                    AlnRecord r;
                    r.query_id  = e.query_id;
                    r.target_id = (int)std::strtol(f1, nullptr, 10);
                    r.raw_score = (int)std::strtol(f2, nullptr, 10);
                    r.fident    = std::strtod(f3, nullptr);
                    r.evalue    = std::strtod(f4, nullptr);
                    r.qstart    = (int)std::strtol(f5, nullptr, 10);
                    r.qend      = (int)std::strtol(f6, nullptr, 10);
                    r.qlen      = (int)std::strtol(f7, nullptr, 10);
                    r.tstart    = (int)std::strtol(f8, nullptr, 10);
                    r.tend      = (int)std::strtol(f9, nullptr, 10);
                    r.tlen      = (int)std::strtol(f10, nullptr, 10);
                    // alnlen = max(|qend-qstart|, |tend-tstart|) + 1  (mmseqs2 Matcher::computeAlnLength)
                    {
                        int dq = r.qend - r.qstart; if (dq < 0) dq = -dq;
                        int dt = r.tend - r.tstart; if (dt < 0) dt = -dt;
                        r.alnlen = (dq > dt ? dq : dt) + 1;
                    }
                    // parse cigar if present: get precise gapopen; mismatch from cigar
                    if (f11 < f11e) {
                        int cigar_alnlen = 0;
                        parse_cigar(f11, f11e, cigar_alnlen, r.gapopen);
                        // use cigar alnlen if available (more accurate)
                        r.alnlen = cigar_alnlen;
                        int nident = (int)(r.fident * r.alnlen + 0.5);
                        r.mismatch = r.alnlen - nident - r.gapopen;
                        if (r.mismatch < 0) r.mismatch = 0;
                    } else {
                        r.gapopen  = 0;
                        // mismatch: mmseqs2 uses min(|qend-qstart|, |tend-tstart|) * (1-seqId) + 0.5
                        // (same as convertalignments.cpp line 447-448, no backtrace branch)
                        int dq = r.qend - r.qstart; if (dq < 0) dq = -dq;
                        int dt = r.tend - r.tstart; if (dt < 0) dt = -dt;
                        float bestMatchEstimate = (float)(dq < dt ? dq : dt);
                        r.mismatch = (int)(bestMatchEstimate * (1.0f - (float)r.fident) + 0.5f);
                        if (r.mismatch < 0) r.mismatch = 0;
                    }

                    const std::string& tname = id2name.count(r.target_id) ? id2name[r.target_id] : std::to_string(r.target_id);
                    write_record(fout, r, qname, tname, fields, sep);
                    total_records++;
                }

                p = (lend < end) ? lend + 1 : end;
            }
        }
    }

    if (!out_path.empty()) std::fclose(fout);

    std::cerr << "Done. Wrote " << total_records << " alignment records.\n";
    return 0;
}
