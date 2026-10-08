#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
namespace fs = std::filesystem;


// ============================================================
// Paths
// ============================================================

//mmseqs相似度结果
const string MMSEQS_FILE =
    "mmseqs_nonself_hits.tsv";

//foldseek相似度结果
const string FOLDSEEK_FILE =
    "foldseek_nonself_hits.tsv";

//是关于蛋白质功能信息
const string ANNOTATION_FILE =
    "Baculoviridae_proteins_name_normalized.csv";

//queries的信息
const string QUERY_FILE =
    "/selected_queries.tsv";

//结构信息的文件    
const string STRUCTURE_MAPPING_FILE =
    "Baculoviridae_structure_mapping.tsv";


//输出的总文件夹
const string OUTPUT_DIR =
    "out/comparison_cpp";


// ============================================================
// Thresholds
// ============================================================

// Weak sequence:
// MMseqs2 identity < 30%
const double SEQ_WEAK_IDENTITY = 0.30;

// Require at least 50% sequence alignment coverage
const double SEQ_MIN_COVERAGE = 0.50;

// Strong structure:
// conservative global TM score
const double STRUCT_TM_THRESHOLD = 0.50;

// Require both sides to have >=50% structural coverage
const double STRUCT_MIN_COVERAGE = 0.50;

// Foldseek homolog probability
const double STRUCT_PROB_THRESHOLD = 0.90;

// Structure quality thresholds
const double PLDDT_THRESHOLD = 70.0;
const double PTM_THRESHOLD   = 0.50;


// ============================================================
// Utility
// ============================================================

double NaN()
{
    return numeric_limits<double>::quiet_NaN();
}


bool isNaN(double x)
{
    return std::isnan(x);
}


string trim(const string& s)
{
    size_t start = s.find_first_not_of(" \t\r\n");

    if (start == string::npos)
        return "";

    size_t end = s.find_last_not_of(" \t\r\n");

    return s.substr(
        start,
        end - start + 1
    );
}


// ============================================================
// TSV splitter
// ============================================================

vector<string> splitTSV(const string& line)
{
    vector<string> fields;

    string field;

    stringstream ss(line);

    while (getline(ss, field, '\t'))
    {
        fields.push_back(field);
    }

    // Preserve final empty field
    if (!line.empty() && line.back() == '\t')
    {
        fields.push_back("");
    }

    return fields;
}


// ============================================================
// CSV parser
//
// Handles quoted CSV such as:
//
// "abc","hello, world","xxx"
//
// ============================================================

vector<string> splitCSV(const string& line)
{
    vector<string> fields;

    string field;

    bool inQuotes = false;

    for (size_t i = 0; i < line.size(); ++i)
    {
        char c = line[i];

        if (c == '"')
        {
            if (
                inQuotes &&
                i + 1 < line.size() &&
                line[i + 1] == '"'
            )
            {
                // Escaped quote: ""
                field += '"';
                ++i;
            }
            else
            {
                inQuotes = !inQuotes;
            }
        }
        else if (
            c == ',' &&
            !inQuotes
        )
        {
            fields.push_back(field);
            field.clear();
        }
        else
        {
            field += c;
        }
    }

    fields.push_back(field);

    return fields;
}


// ============================================================
// Header helper
// ============================================================

unordered_map<string, size_t>
buildHeaderMap(
    const vector<string>& headers
)
{
    unordered_map<string, size_t> result;

    for (size_t i = 0; i < headers.size(); ++i)
    {
        string name = trim(headers[i]);

        // Remove UTF-8 BOM if present
        if (
            i == 0 &&
            name.size() >= 3 &&
            (unsigned char)name[0] == 0xEF &&
            (unsigned char)name[1] == 0xBB &&
            (unsigned char)name[2] == 0xBF
        )
        {
            name = name.substr(3);
        }

        result[name] = i;
    }

    return result;
}


string getField(
    const vector<string>& fields,
    const unordered_map<string, size_t>& header,
    const string& name
)
{
    auto it = header.find(name);

    if (it == header.end())
        return "";

    size_t idx = it->second;

    if (idx >= fields.size())
        return "";

    return trim(fields[idx]);
}


double toDouble(
    const string& s
)
{
    string x = trim(s);

    if (x.empty())
        return NaN();

    try
    {
        return stod(x);
    }
    catch (...)
    {
        return NaN();
    }
}


// ============================================================
// Pair key
// ============================================================

string pairKey(
    const string& query,
    const string& target
)
{
    return query + "\t" + target;
}


// ============================================================
// Sequence result
// ============================================================

struct SequenceHit
{
    bool exists = false;

    string query;
    string target;

    double fident   = NaN();
    double alnlen   = NaN();
    double mismatch = NaN();
    double gapopen  = NaN();

    double qstart = NaN();
    double qend   = NaN();

    double tstart = NaN();
    double tend   = NaN();

    double evalue = NaN();
    double bits   = NaN();

    double qlen = NaN();
    double tlen = NaN();

    double qcov = NaN();
    double tcov = NaN();
};


// ============================================================
// Foldseek result
// ============================================================

struct StructureHit
{
    bool exists = false;

    string query;
    string target;

    double fident   = NaN();
    double alnlen   = NaN();
    double mismatch = NaN();
    double gapopen  = NaN();

    double qstart = NaN();
    double qend   = NaN();

    double tstart = NaN();
    double tend   = NaN();

    double evalue = NaN();
    double bits   = NaN();

    double qlen = NaN();
    double tlen = NaN();

    double qcov = NaN();
    double tcov = NaN();

    double alntmscore = NaN();
    double qtmscore   = NaN();
    double ttmscore   = NaN();

    double lddt = NaN();
    double prob = NaN();
};


// ============================================================
// Protein annotation
// ============================================================

struct Annotation
{
    string id;

    string species;

    string viro3dName;

    string standardizedName;

    string genbankProteinID;

    string uniprotID;

    double length = NaN();

    double pLDDT = NaN();

    double pTM = NaN();
};


// ============================================================
// Query metadata
// ============================================================

struct QueryInfo
{
    string standardizedFunction;

    int groupRank = -1;
};


// ============================================================
// Structure mapping
// ============================================================

struct StructureInfo
{
    bool available = false;

    string status;

    string pdbPath;

    double pLDDT = NaN();

    double pTM = NaN();
};


// ============================================================
// Final pair record
// ============================================================

struct PairRecord
{
    string query;
    string target;

    SequenceHit seq;

    StructureHit str;

    Annotation queryAnno;
    Annotation targetAnno;

    QueryInfo queryInfo;

    StructureInfo targetStructure;

    bool weakSequence = false;

    bool strongStructure = false;

    bool queryQualityOK = false;

    bool targetQualityOK = false;

    bool bothQualityOK = false;

    bool weakSeqStrongStruct = false;

    bool seqUndetectedStrongStruct = false;

    string category;

    double seqMinCov = NaN();

    double structMinCov = NaN();

    double structGlobalTM = NaN();

    double candidateScore = NaN();
};


// ============================================================
// Read MMseqs2
// ============================================================

unordered_map<string, SequenceHit>
readMMseqs(
    const string& path
)
{
    ifstream fin(path);

    if (!fin)
    {
        throw runtime_error(
            "Cannot open MMseqs2 file: " + path
        );
    }

    string line;

    getline(fin, line);

    auto headers = splitTSV(line);

    auto header = buildHeaderMap(headers);

    unordered_map<string, SequenceHit> result;

    size_t n = 0;

    while (getline(fin, line))
    {
        if (trim(line).empty())
            continue;

        auto fields = splitTSV(line);

        SequenceHit h;

        h.exists = true;

        h.query = getField(
            fields,
            header,
            "query"
        );

        h.target = getField(
            fields,
            header,
            "target"
        );

        h.fident = toDouble(
            getField(fields, header, "fident")
        );

        h.alnlen = toDouble(
            getField(fields, header, "alnlen")
        );

        h.mismatch = toDouble(
            getField(fields, header, "mismatch")
        );

        h.gapopen = toDouble(
            getField(fields, header, "gapopen")
        );

        h.qstart = toDouble(
            getField(fields, header, "qstart")
        );

        h.qend = toDouble(
            getField(fields, header, "qend")
        );

        h.tstart = toDouble(
            getField(fields, header, "tstart")
        );

        h.tend = toDouble(
            getField(fields, header, "tend")
        );

        h.evalue = toDouble(
            getField(fields, header, "evalue")
        );

        h.bits = toDouble(
            getField(fields, header, "bits")
        );

        h.qlen = toDouble(
            getField(fields, header, "qlen")
        );

        h.tlen = toDouble(
            getField(fields, header, "tlen")
        );

        h.qcov = toDouble(
            getField(fields, header, "qcov")
        );

        h.tcov = toDouble(
            getField(fields, header, "tcov")
        );

        string key = pairKey(
            h.query,
            h.target
        );

        auto it = result.find(key);

        // If duplicate pair exists:
        // keep lower E-value;
        // if equal, keep higher bit score.
        if (it == result.end())
        {
            result[key] = h;
        }
        else
        {
            SequenceHit& old = it->second;

            bool replace = false;

            if (
                !isNaN(h.evalue) &&
                (
                    isNaN(old.evalue) ||
                    h.evalue < old.evalue
                )
            )
            {
                replace = true;
            }
            else if (
                h.evalue == old.evalue &&
                !isNaN(h.bits) &&
                (
                    isNaN(old.bits) ||
                    h.bits > old.bits
                )
            )
            {
                replace = true;
            }

            if (replace)
                old = h;
        }

        ++n;
    }

    cout
        << "MMseqs2 rows read: "
        << n
        << "\n";

    cout
        << "MMseqs2 unique pairs: "
        << result.size()
        << "\n";

    return result;
}


// ============================================================
// Read Foldseek
// ============================================================

unordered_map<string, StructureHit>
readFoldseek(
    const string& path
)
{
    ifstream fin(path);

    if (!fin)
    {
        throw runtime_error(
            "Cannot open Foldseek file: " + path
        );
    }

    string line;

    getline(fin, line);

    auto headers = splitTSV(line);

    auto header = buildHeaderMap(headers);

    unordered_map<string, StructureHit> result;

    size_t n = 0;

    while (getline(fin, line))
    {
        if (trim(line).empty())
            continue;

        auto fields = splitTSV(line);

        StructureHit h;

        h.exists = true;

        h.query = getField(
            fields,
            header,
            "query"
        );

        h.target = getField(
            fields,
            header,
            "target"
        );

        h.fident = toDouble(
            getField(fields, header, "fident")
        );

        h.alnlen = toDouble(
            getField(fields, header, "alnlen")
        );

        h.mismatch = toDouble(
            getField(fields, header, "mismatch")
        );

        h.gapopen = toDouble(
            getField(fields, header, "gapopen")
        );

        h.qstart = toDouble(
            getField(fields, header, "qstart")
        );

        h.qend = toDouble(
            getField(fields, header, "qend")
        );

        h.tstart = toDouble(
            getField(fields, header, "tstart")
        );

        h.tend = toDouble(
            getField(fields, header, "tend")
        );

        h.evalue = toDouble(
            getField(fields, header, "evalue")
        );

        h.bits = toDouble(
            getField(fields, header, "bits")
        );

        h.qlen = toDouble(
            getField(fields, header, "qlen")
        );

        h.tlen = toDouble(
            getField(fields, header, "tlen")
        );

        h.qcov = toDouble(
            getField(fields, header, "qcov")
        );

        h.tcov = toDouble(
            getField(fields, header, "tcov")
        );

        h.alntmscore = toDouble(
            getField(fields, header, "alntmscore")
        );

        h.qtmscore = toDouble(
            getField(fields, header, "qtmscore")
        );

        h.ttmscore = toDouble(
            getField(fields, header, "ttmscore")
        );

        h.lddt = toDouble(
            getField(fields, header, "lddt")
        );

        h.prob = toDouble(
            getField(fields, header, "prob")
        );

        string key = pairKey(
            h.query,
            h.target
        );

        auto it = result.find(key);

        if (it == result.end())
        {
            result[key] = h;
        }
        else
        {
            StructureHit& old = it->second;

            double newTM = NaN();
            double oldTM = NaN();

            if (
                !isNaN(h.qtmscore) &&
                !isNaN(h.ttmscore)
            )
            {
                newTM = min(
                    h.qtmscore,
                    h.ttmscore
                );
            }

            if (
                !isNaN(old.qtmscore) &&
                !isNaN(old.ttmscore)
            )
            {
                oldTM = min(
                    old.qtmscore,
                    old.ttmscore
                );
            }

            bool replace = false;

            if (
                !isNaN(newTM) &&
                (
                    isNaN(oldTM) ||
                    newTM > oldTM
                )
            )
            {
                replace = true;
            }
            else if (
                newTM == oldTM &&
                h.prob > old.prob
            )
            {
                replace = true;
            }

            if (replace)
                old = h;
        }

        ++n;
    }

    cout
        << "Foldseek rows read: "
        << n
        << "\n";

    cout
        << "Foldseek unique pairs: "
        << result.size()
        << "\n";

    return result;
}


// ============================================================
// Read annotation CSV
// ============================================================

unordered_map<string, Annotation>
readAnnotations(
    const string& path
)
{
    ifstream fin(path);

    if (!fin)
    {
        throw runtime_error(
            "Cannot open annotation CSV: " + path
        );
    }

    string line;

    getline(fin, line);

    auto headers = splitCSV(line);

    auto header = buildHeaderMap(headers);

    unordered_map<string, Annotation> result;

    size_t n = 0;

    while (getline(fin, line))
    {
        if (trim(line).empty())
            continue;

        auto fields = splitCSV(line);

        Annotation a;

        a.id = getField(
            fields,
            header,
            "Viro3D ID"
        );

        if (a.id.empty())
            continue;

        a.species = getField(
            fields,
            header,
            "ICTV Species"
        );

        a.viro3dName = getField(
            fields,
            header,
            "Viro3D Name"
        );

        a.standardizedName = getField(
            fields,
            header,
            "Standardized Viro3D Name"
        );

        a.genbankProteinID = getField(
            fields,
            header,
            "GenBank Protein ID"
        );

        a.uniprotID = getField(
            fields,
            header,
            "UniProt ID"
        );

        a.length = toDouble(
            getField(
                fields,
                header,
                "Protein Length"
            )
        );

        a.pLDDT = toDouble(
            getField(
                fields,
                header,
                "ColabFold pLDDT"
            )
        );

        a.pTM = toDouble(
            getField(
                fields,
                header,
                "ColabFold pTM"
            )
        );

        if (!result.count(a.id))
        {
            result[a.id] = a;
        }

        ++n;
    }

    cout
        << "Annotation rows read: "
        << n
        << "\n";

    cout
        << "Unique annotated proteins: "
        << result.size()
        << "\n";

    return result;
}


// ============================================================
// Read query metadata
// ============================================================

unordered_map<string, QueryInfo>
readQueryInfo(
    const string& path
)
{
    ifstream fin(path);

    if (!fin)
    {
        throw runtime_error(
            "Cannot open query TSV: " + path
        );
    }

    string line;

    getline(fin, line);

    auto headers = splitTSV(line);

    auto header = buildHeaderMap(headers);

    unordered_map<string, QueryInfo> result;

    while (getline(fin, line))
    {
        if (trim(line).empty())
            continue;

        auto fields = splitTSV(line);

        string id = getField(
            fields,
            header,
            "Viro3D_ID"
        );

        if (id.empty())
            continue;

        QueryInfo q;

        q.standardizedFunction = getField(
            fields,
            header,
            "Standardized_Function"
        );

        string rank = getField(
            fields,
            header,
            "Query_Group_Rank"
        );

        if (!rank.empty())
        {
            try
            {
                q.groupRank = stoi(rank);
            }
            catch (...)
            {
                q.groupRank = -1;
            }
        }

        result[id] = q;
    }

    cout
        << "Selected queries: "
        << result.size()
        << "\n";

    return result;
}


// ============================================================
// Read structure mapping
// ============================================================

unordered_map<string, StructureInfo>
readStructureMapping(
    const string& path
)
{
    ifstream fin(path);

    if (!fin)
    {
        throw runtime_error(
            "Cannot open structure mapping: " + path
        );
    }

    string line;

    getline(fin, line);

    auto headers = splitTSV(line);

    auto header = buildHeaderMap(headers);

    unordered_map<string, StructureInfo> result;

    while (getline(fin, line))
    {
        if (trim(line).empty())
            continue;

        auto fields = splitTSV(line);

        string id = getField(
            fields,
            header,
            "Viro3D_ID"
        );

        if (id.empty())
            continue;

        StructureInfo s;

        s.status = getField(
            fields,
            header,
            "Structure_Status"
        );

        s.pdbPath = getField(
            fields,
            header,
            "Dataset_PDB_Path"
        );

        s.pLDDT = toDouble(
            getField(
                fields,
                header,
                "ColabFold_pLDDT"
            )
        );

        s.pTM = toDouble(
            getField(
                fields,
                header,
                "ColabFold_pTM"
            )
        );

        // Structure is usable only when the path exists.
        s.available =
            !s.pdbPath.empty() &&
            fs::exists(s.pdbPath);

        result[id] = s;
    }

    cout
        << "Structure mapping proteins: "
        << result.size()
        << "\n";

    return result;
}


// ============================================================
// Safe output helper
// ============================================================

string num(double x)
{
    if (isNaN(x))
        return "NA";

    ostringstream oss;

    oss
        << setprecision(8)
        << x;

    return oss.str();
}


string boolStr(bool x)
{
    return x ? "1" : "0";
}


// ============================================================
// Build final records
// ============================================================

vector<PairRecord>
buildRecords(
    const unordered_map<string, SequenceHit>& seqHits,
    const unordered_map<string, StructureHit>& structHits,
    const unordered_map<string, Annotation>& annotations,
    const unordered_map<string, QueryInfo>& queryInfo,
    const unordered_map<string, StructureInfo>& structures
)
{
    unordered_set<string> allKeys;

    for (const auto& kv : seqHits)
    {
        allKeys.insert(kv.first);
    }

    for (const auto& kv : structHits)
    {
        allKeys.insert(kv.first);
    }

    vector<PairRecord> records;

    records.reserve(
        allKeys.size()
    );

    for (const string& key : allKeys)
    {
        PairRecord r;

        auto sIt = seqHits.find(key);

        if (sIt != seqHits.end())
        {
            r.seq = sIt->second;

            r.query = r.seq.query;

            r.target = r.seq.target;
        }

        auto fIt = structHits.find(key);

        if (fIt != structHits.end())
        {
            r.str = fIt->second;

            if (r.query.empty())
            {
                r.query = r.str.query;

                r.target = r.str.target;
            }
        }

        // ----------------------------------------------------
        // Annotation
        // ----------------------------------------------------

        auto qa = annotations.find(
            r.query
        );

        if (qa != annotations.end())
        {
            r.queryAnno = qa->second;
        }

        auto ta = annotations.find(
            r.target
        );

        if (ta != annotations.end())
        {
            r.targetAnno = ta->second;
        }

        // ----------------------------------------------------
        // Query metadata
        // ----------------------------------------------------

        auto qi = queryInfo.find(
            r.query
        );

        if (qi != queryInfo.end())
        {
            r.queryInfo = qi->second;
        }

        // ----------------------------------------------------
        // Target structure availability
        // ----------------------------------------------------

        auto st = structures.find(
            r.target
        );

        if (st != structures.end())
        {
            r.targetStructure = st->second;
        }

        // ----------------------------------------------------
        // Sequence minimum coverage
        // ----------------------------------------------------

        if (
            r.seq.exists &&
            !isNaN(r.seq.qcov) &&
            !isNaN(r.seq.tcov)
        )
        {
            r.seqMinCov = min(
                r.seq.qcov,
                r.seq.tcov
            );
        }

        // ----------------------------------------------------
        // Structure minimum coverage
        // ----------------------------------------------------

        if (
            r.str.exists &&
            !isNaN(r.str.qcov) &&
            !isNaN(r.str.tcov)
        )
        {
            r.structMinCov = min(
                r.str.qcov,
                r.str.tcov
            );
        }

        // ----------------------------------------------------
        // Conservative global TM
        //
        // Use min(qTM, tTM)
        // ----------------------------------------------------

        if (
            r.str.exists &&
            !isNaN(r.str.qtmscore) &&
            !isNaN(r.str.ttmscore)
        )
        {
            r.structGlobalTM = min(
                r.str.qtmscore,
                r.str.ttmscore
            );
        }

        // ----------------------------------------------------
        // Weak sequence
        // ----------------------------------------------------

        r.weakSequence =
            r.seq.exists
            &&
            !isNaN(r.seq.fident)
            &&
            !isNaN(r.seqMinCov)
            &&
            r.seq.fident < SEQ_WEAK_IDENTITY
            &&
            r.seqMinCov >= SEQ_MIN_COVERAGE;

        // ----------------------------------------------------
        // Strong structure
        // ----------------------------------------------------

        r.strongStructure =
            r.str.exists
            &&
            !isNaN(r.structGlobalTM)
            &&
            !isNaN(r.structMinCov)
            &&
            !isNaN(r.str.prob)
            &&
            r.structGlobalTM >= STRUCT_TM_THRESHOLD
            &&
            r.structMinCov >= STRUCT_MIN_COVERAGE
            &&
            r.str.prob >= STRUCT_PROB_THRESHOLD;

        // ----------------------------------------------------
        // Structure quality
        // ----------------------------------------------------

        r.queryQualityOK =
            !isNaN(r.queryAnno.pLDDT)
            &&
            !isNaN(r.queryAnno.pTM)
            &&
            r.queryAnno.pLDDT >= PLDDT_THRESHOLD
            &&
            r.queryAnno.pTM >= PTM_THRESHOLD;

        r.targetQualityOK =
            !isNaN(r.targetAnno.pLDDT)
            &&
            !isNaN(r.targetAnno.pTM)
            &&
            r.targetAnno.pLDDT >= PLDDT_THRESHOLD
            &&
            r.targetAnno.pTM >= PTM_THRESHOLD;

        r.bothQualityOK =
            r.queryQualityOK
            &&
            r.targetQualityOK;

        // ----------------------------------------------------
        // Main category
        // ----------------------------------------------------

        if (
            r.seq.exists &&
            r.str.exists
        )
        {
            r.category = "Both";
        }
        else if (
            r.seq.exists &&
            !r.str.exists
        )
        {
            if (
                r.targetStructure.available
            )
            {
                r.category =
                    "Sequence_only";
            }
            else
            {
                r.category =
                    "Sequence_hit_structure_unavailable";
            }
        }
        else if (
            !r.seq.exists &&
            r.str.exists
        )
        {
            r.category =
                "Structure_only";
        }
        else
        {
            r.category = "Other";
        }

        // ----------------------------------------------------
        // Special categories
        // ----------------------------------------------------

        r.weakSeqStrongStruct =
            r.weakSequence
            &&
            r.strongStructure;

        r.seqUndetectedStrongStruct =
            !r.seq.exists
            &&
            r.strongStructure;

        // ----------------------------------------------------
        // Candidate ranking score
        //
        // This is only for ranking.
        // It is NOT a biological similarity metric.
        // ----------------------------------------------------

        if (r.strongStructure)
        {
            double seqWeakness =
                r.seq.exists &&
                !isNaN(r.seq.fident)
                ?
                (1.0 - r.seq.fident)
                :
                1.0;

            double score = 0.0;

            score +=
                3.0 *
                (
                    isNaN(r.structGlobalTM)
                    ?
                    0.0
                    :
                    r.structGlobalTM
                );

            score +=
                2.0 *
                (
                    isNaN(r.structMinCov)
                    ?
                    0.0
                    :
                    r.structMinCov
                );

            score +=
                (
                    isNaN(r.str.prob)
                    ?
                    0.0
                    :
                    r.str.prob
                );

            score += seqWeakness;

            if (r.bothQualityOK)
            {
                score += 1.0;
            }

            r.candidateScore = score;
        }

        records.push_back(
            r
        );
    }

    return records;
}


// ============================================================
// Output header
// ============================================================

void writeHeader(
    ofstream& out
)
{
    out
        << "query\t"
        << "target\t"

        << "query_species\t"
        << "target_species\t"

        << "query_function\t"
        << "query_annotation\t"
        << "target_annotation\t"

        << "seq_hit\t"
        << "seq_fident\t"
        << "seq_identity_percent\t"
        << "seq_evalue\t"
        << "seq_bits\t"
        << "seq_qcov\t"
        << "seq_tcov\t"
        << "seq_min_cov\t"

        << "struct_hit\t"
        << "struct_fident\t"
        << "struct_evalue\t"
        << "struct_bits\t"
        << "struct_qcov\t"
        << "struct_tcov\t"
        << "struct_min_cov\t"

        << "struct_alntmscore\t"
        << "struct_qtmscore\t"
        << "struct_ttmscore\t"
        << "struct_global_tm\t"
        << "struct_lddt\t"
        << "struct_prob\t"

        << "query_pLDDT\t"
        << "query_pTM\t"
        << "target_pLDDT\t"
        << "target_pTM\t"

        << "target_structure_available\t"
        << "target_structure_status\t"
        << "target_pdb_path\t"

        << "hit_category\t"

        << "weak_sequence\t"
        << "strong_structure\t"

        << "weak_sequence_strong_structure\t"
        << "sequence_undetected_strong_structure\t"

        << "both_structure_quality_ok\t"

        << "candidate_score"

        << "\n";
}


// ============================================================
// Output one record
// ============================================================

void writeRecord(
    ofstream& out,
    const PairRecord& r
)
{
    double identityPercent =
        (
            r.seq.exists &&
            !isNaN(r.seq.fident)
        )
        ?
        r.seq.fident * 100.0
        :
        NaN();

    out
        << r.query << "\t"
        << r.target << "\t"

        << r.queryAnno.species << "\t"
        << r.targetAnno.species << "\t"

        << r.queryInfo.standardizedFunction << "\t"
        << r.queryAnno.standardizedName << "\t"
        << r.targetAnno.standardizedName << "\t"

        << boolStr(r.seq.exists) << "\t"
        << num(r.seq.fident) << "\t"
        << num(identityPercent) << "\t"
        << num(r.seq.evalue) << "\t"
        << num(r.seq.bits) << "\t"
        << num(r.seq.qcov) << "\t"
        << num(r.seq.tcov) << "\t"
        << num(r.seqMinCov) << "\t"

        << boolStr(r.str.exists) << "\t"
        << num(r.str.fident) << "\t"
        << num(r.str.evalue) << "\t"
        << num(r.str.bits) << "\t"
        << num(r.str.qcov) << "\t"
        << num(r.str.tcov) << "\t"
        << num(r.structMinCov) << "\t"

        << num(r.str.alntmscore) << "\t"
        << num(r.str.qtmscore) << "\t"
        << num(r.str.ttmscore) << "\t"
        << num(r.structGlobalTM) << "\t"
        << num(r.str.lddt) << "\t"
        << num(r.str.prob) << "\t"

        << num(r.queryAnno.pLDDT) << "\t"
        << num(r.queryAnno.pTM) << "\t"
        << num(r.targetAnno.pLDDT) << "\t"
        << num(r.targetAnno.pTM) << "\t"

        << boolStr(
            r.targetStructure.available
        )
        << "\t"

        << r.targetStructure.status
        << "\t"

        << r.targetStructure.pdbPath
        << "\t"

        << r.category << "\t"

        << boolStr(
            r.weakSequence
        )
        << "\t"

        << boolStr(
            r.strongStructure
        )
        << "\t"

        << boolStr(
            r.weakSeqStrongStruct
        )
        << "\t"

        << boolStr(
            r.seqUndetectedStrongStruct
        )
        << "\t"

        << boolStr(
            r.bothQualityOK
        )
        << "\t"

        << num(
            r.candidateScore
        )

        << "\n";
}


// ============================================================
// Save a subset
// ============================================================

template <typename Predicate>
void writeSubset(
    const string& filename,
    const vector<PairRecord>& records,
    Predicate keep
)
{
    ofstream out(filename);

    if (!out)
    {
        throw runtime_error(
            "Cannot write: " + filename
        );
    }

    writeHeader(out);

    for (const auto& r : records)
    {
        if (keep(r))
        {
            writeRecord(
                out,
                r
            );
        }
    }
}


// ============================================================
// Main
// ============================================================

int main()
{
    try
    {
        cout
            << "============================================================\n"
            << "MMseqs2 + Foldseek integration\n"
            << "============================================================\n";

        fs::create_directories(
            OUTPUT_DIR
        );

        // ----------------------------------------------------
        // 1. Read input files
        // ----------------------------------------------------

        cout << "\n[1] Reading MMseqs2\n";

        auto seqHits =
            readMMseqs(
                MMSEQS_FILE
            );

        cout << "\n[2] Reading Foldseek\n";

        auto structHits =
            readFoldseek(
                FOLDSEEK_FILE
            );

        cout << "\n[3] Reading annotation\n";

        auto annotations =
            readAnnotations(
                ANNOTATION_FILE
            );

        cout << "\n[4] Reading query metadata\n";

        auto queryInfo =
            readQueryInfo(
                QUERY_FILE
            );

        cout << "\n[5] Reading structure mapping\n";

        auto structures =
            readStructureMapping(
                STRUCTURE_MAPPING_FILE
            );

        // ----------------------------------------------------
        // 2. Integrate
        // ----------------------------------------------------

        cout
            << "\n[6] Integrating query-target pairs\n";

        vector<PairRecord> records =
            buildRecords(
                seqHits,
                structHits,
                annotations,
                queryInfo,
                structures
            );

        cout
            << "Integrated unique pairs: "
            << records.size()
            << "\n";

        // ----------------------------------------------------
        // 3. Sort all pairs
        //
        // query first, then target
        // ----------------------------------------------------

        sort(
            records.begin(),
            records.end(),
            [](
                const PairRecord& a,
                const PairRecord& b
            )
            {
                if (a.query != b.query)
                    return a.query < b.query;

                return a.target < b.target;
            }
        );

        // ----------------------------------------------------
        // 4. Save integrated table
        // ----------------------------------------------------

        string integratedFile =
            OUTPUT_DIR
            + "/integrated_all_pairs.tsv";

        {
            ofstream out(
                integratedFile
            );

            writeHeader(out);

            for (const auto& r : records)
            {
                writeRecord(
                    out,
                    r
                );
            }
        }

        // ----------------------------------------------------
        // 5. Save categories
        // ----------------------------------------------------

        writeSubset(
            OUTPUT_DIR
            + "/both_hits.tsv",
            records,
            [](
                const PairRecord& r
            )
            {
                return
                    r.category == "Both";
            }
        );

        writeSubset(
            OUTPUT_DIR
            + "/sequence_only.tsv",
            records,
            [](
                const PairRecord& r
            )
            {
                return
                    r.category
                    == "Sequence_only";
            }
        );

        writeSubset(
            OUTPUT_DIR
            + "/structure_only.tsv",
            records,
            [](
                const PairRecord& r
            )
            {
                return
                    r.category
                    == "Structure_only";
            }
        );

        writeSubset(
            OUTPUT_DIR
            + "/sequence_hit_structure_unavailable.tsv",
            records,
            [](
                const PairRecord& r
            )
            {
                return
                    r.category
                    ==
                    "Sequence_hit_structure_unavailable";
            }
        );

        writeSubset(
            OUTPUT_DIR
            + "/weak_sequence_strong_structure.tsv",
            records,
            [](
                const PairRecord& r
            )
            {
                return
                    r.weakSeqStrongStruct;
            }
        );

        writeSubset(
            OUTPUT_DIR
            + "/sequence_undetected_strong_structure.tsv",
            records,
            [](
                const PairRecord& r
            )
            {
                return
                    r.seqUndetectedStrongStruct;
            }
        );

        // ----------------------------------------------------
        // 6. Representative candidates
        //
        // weak-sequence + strong-structure
        // OR
        // no-MMseqs2 + strong-structure
        //
        // AND
        // both predicted structures are reliable
        // ----------------------------------------------------

        vector<PairRecord> candidates;

        for (const auto& r : records)
        {
            bool interesting =
                r.weakSeqStrongStruct
                ||
                r.seqUndetectedStrongStruct;

            if (
                interesting
                &&
                r.bothQualityOK
            )
            {
                candidates.push_back(
                    r
                );
            }
        }

        sort(
            candidates.begin(),
            candidates.end(),
            [](
                const PairRecord& a,
                const PairRecord& b
            )
            {
                double sa =
                    isNaN(a.candidateScore)
                    ?
                    -1.0
                    :
                    a.candidateScore;

                double sb =
                    isNaN(b.candidateScore)
                    ?
                    -1.0
                    :
                    b.candidateScore;

                return sa > sb;
            }
        );

        string candidateFile =
            OUTPUT_DIR
            + "/representative_candidates.tsv";

        {
            ofstream out(
                candidateFile
            );

            writeHeader(out);

            for (const auto& r : candidates)
            {
                writeRecord(
                    out,
                    r
                );
            }
        }

        // ----------------------------------------------------
        // 7. Statistics
        // ----------------------------------------------------

        map<string, size_t> categoryCount;

        size_t weakStrongCount = 0;

        size_t undetectedStrongCount = 0;

        size_t qualityCandidateCount = 0;

        unordered_set<string> seqQueries;

        unordered_set<string> structQueries;

        for (const auto& r : records)
        {
            categoryCount[
                r.category
            ]++;

            if (
                r.weakSeqStrongStruct
            )
            {
                ++weakStrongCount;
            }

            if (
                r.seqUndetectedStrongStruct
            )
            {
                ++undetectedStrongCount;
            }

            if (
                (
                    r.weakSeqStrongStruct
                    ||
                    r.seqUndetectedStrongStruct
                )
                &&
                r.bothQualityOK
            )
            {
                ++qualityCandidateCount;
            }

            if (r.seq.exists)
            {
                seqQueries.insert(
                    r.query
                );
            }

            if (r.str.exists)
            {
                structQueries.insert(
                    r.query
                );
            }
        }

        // ----------------------------------------------------
        // Structure availability statistics
        // ----------------------------------------------------

        size_t totalStructures = 0;

        size_t availableStructures = 0;

        size_t unavailableStructures = 0;

        for (
            const auto& kv :
            structures
        )
        {
            ++totalStructures;

            if (
                kv.second.available
            )
            {
                ++availableStructures;
            }
            else
            {
                ++unavailableStructures;
            }
        }

        // ----------------------------------------------------
        // Write summary
        // ----------------------------------------------------

        string summaryFile =
            OUTPUT_DIR
            + "/summary_statistics.txt";

        ofstream summary(
            summaryFile
        );

        summary
            << "Baculoviridae MMseqs2 vs Foldseek comparison\n"
            << "============================================================\n\n";

        summary
            << "MMseqs2 unique non-self pairs: "
            << seqHits.size()
            << "\n";

        summary
            << "Foldseek unique non-self pairs: "
            << structHits.size()
            << "\n";

        summary
            << "Integrated unique pairs: "
            << records.size()
            << "\n\n";

        summary
            << "Queries with MMseqs2 hits: "
            << seqQueries.size()
            << "\n";

        summary
            << "Queries with Foldseek hits: "
            << structQueries.size()
            << "\n\n";

        summary
            << "Hit categories:\n";

        for (
            const auto& kv :
            categoryCount
        )
        {
            summary
                << "  "
                << kv.first
                << ": "
                << kv.second
                << "\n";
        }

        summary << "\n";

        summary
            << "Weak sequence + strong structure: "
            << weakStrongCount
            << "\n";

        summary
            << "MMseqs2-undetected + strong structure: "
            << undetectedStrongCount
            << "\n";

        summary
            << "High-quality representative candidates: "
            << qualityCandidateCount
            << "\n\n";

        summary
            << "Structure availability:\n";

        summary
            << "  Mapping records: "
            << totalStructures
            << "\n";

        summary
            << "  Usable structure files: "
            << availableStructures
            << "\n";

        summary
            << "  Missing/unavailable structure files: "
            << unavailableStructures
            << "\n\n";

        summary
            << "Thresholds:\n";

        summary
            << "  Weak sequence identity < "
            << SEQ_WEAK_IDENTITY
            << "\n";

        summary
            << "  Minimum sequence coverage >= "
            << SEQ_MIN_COVERAGE
            << "\n";

        summary
            << "  Strong structure min(qTM,tTM) >= "
            << STRUCT_TM_THRESHOLD
            << "\n";

        summary
            << "  Minimum structure coverage >= "
            << STRUCT_MIN_COVERAGE
            << "\n";

        summary
            << "  Foldseek probability >= "
            << STRUCT_PROB_THRESHOLD
            << "\n";

        summary
            << "  pLDDT >= "
            << PLDDT_THRESHOLD
            << "\n";

        summary
            << "  pTM >= "
            << PTM_THRESHOLD
            << "\n";

        summary.close();

        // ----------------------------------------------------
        // Console summary
        // ----------------------------------------------------

        cout
            << "\n============================================================\n"
            << "Analysis completed\n"
            << "============================================================\n";

        cout
            << "MMseqs2 unique hits       : "
            << seqHits.size()
            << "\n";

        cout
            << "Foldseek unique hits      : "
            << structHits.size()
            << "\n";

        cout
            << "Integrated unique pairs   : "
            << records.size()
            << "\n\n";

        for (
            const auto& kv :
            categoryCount
        )
        {
            cout
                << left
                << setw(38)
                << kv.first
                << ": "
                << kv.second
                << "\n";
        }

        cout
            << "\nWeak sequence + strong structure : "
            << weakStrongCount
            << "\n";

        cout
            << "No MMseqs2 + strong structure     : "
            << undetectedStrongCount
            << "\n";

        cout
            << "High-quality candidates           : "
            << qualityCandidateCount
            << "\n";

        cout
            << "\nStructure mapping records         : "
            << totalStructures
            << "\n";

        cout
            << "Usable structures                 : "
            << availableStructures
            << "\n";

        cout
            << "Unavailable structures            : "
            << unavailableStructures
            << "\n";

        cout
            << "\nOutput directory:\n"
            << OUTPUT_DIR
            << "\n";

        cout
            << "\nMost important outputs:\n";

        cout
            << integratedFile
            << "\n";

        cout
            << OUTPUT_DIR
            << "/weak_sequence_strong_structure.tsv\n";

        cout
            << OUTPUT_DIR
            << "/sequence_undetected_strong_structure.tsv\n";

        cout
            << candidateFile
            << "\n";

        cout
            << summaryFile
            << "\n";

        cout
            << "============================================================\n";

        return 0;
    }
    catch (
        const exception& e
    )
    {
        cerr
            << "ERROR: "
            << e.what()
            << "\n";

        return 1;
    }
}