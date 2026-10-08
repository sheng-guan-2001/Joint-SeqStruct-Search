#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
namespace fs = std::filesystem;


// ============================================================
// 去掉两端空格
// ============================================================
string trim(const string& s) {

    size_t start = 0;

    while (start < s.size() &&
           isspace(static_cast<unsigned char>(s[start]))) {
        ++start;
    }

    size_t end = s.size();

    while (end > start &&
           isspace(static_cast<unsigned char>(s[end - 1]))) {
        --end;
    }

    return s.substr(start, end - start);
}


// ============================================================
// CSV parser
// 支持带逗号和引号的CSV字段
// ============================================================
vector<string> parseCSVLine(const string& line) {

    vector<string> fields;
    string current;

    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {

        char c = line[i];

        if (c == '"') {

            if (in_quotes &&
                i + 1 < line.size() &&
                line[i + 1] == '"') {

                current += '"';
                ++i;

            } else {

                in_quotes = !in_quotes;
            }

        } else if (c == ',' && !in_quotes) {

            fields.push_back(current);
            current.clear();

        } else {

            current += c;
        }
    }

    fields.push_back(current);

    return fields;
}


// ============================================================
// 判断字符串是否以某个suffix结尾
// ============================================================
bool endsWith(
    const string& str,
    const string& suffix
) {

    if (str.size() < suffix.size()) {
        return false;
    }

    return str.compare(
        str.size() - suffix.size(),
        suffix.size(),
        suffix
    ) == 0;
}


// ============================================================
// 从PDB文件名获得Viro3D ID
//
// 例如：
// AAA02739.1.1.10_7195_relaxed.pdb
//
// ->
//
// AAA02739.1.1.10_7195
// ============================================================
string pdbFilenameToViro3DID(
    const string& filename
) {

    const string suffix =
        "_relaxed.pdb";

    if (!endsWith(filename, suffix)) {
        return "";
    }

    return filename.substr(
        0,
        filename.size() - suffix.size()
    );
}


// ============================================================
// 创建软链接
// ============================================================
bool createSymlinkSafe(
    const fs::path& source,
    const fs::path& destination
) {

    try {

        // 如果已经存在
        if (
            fs::exists(destination) ||
            fs::is_symlink(destination)
        ) {

            fs::remove(destination);
        }

        fs::create_symlink(
            source,
            destination
        );

        return true;

    } catch (const fs::filesystem_error& e) {

        cerr
            << "WARNING: cannot create symlink:\n"
            << destination
            << "\nReason: "
            << e.what()
            << endl;

        return false;
    }
}


// ============================================================
// 主函数
// ============================================================
int main() {

    // ========================================================
    // 输入：Baculoviridae protein metadata
    // ========================================================

    const string protein_csv =
        "Baculoviridae_proteins.csv";


    // ========================================================
    // 输入：所有Viro3D ColabFold PDB
    // ========================================================

    const string pdb_root =
        "raw_unzip/colabfold_pdb/colabfold_pdb";


    // ========================================================
    // 输出目录
    // ========================================================

    const string output_root =
        "processed/Baculoviridae";


    const string structure_dir =
        output_root +
        "/structures";


    fs::create_directories(
        structure_dir
    );


    // ========================================================
    // 输出文件
    // ========================================================

    const string mapping_file =
        output_root +
        "/Baculoviridae_structure_mapping.tsv";


    const string missing_file =
        output_root +
        "/Baculoviridae_missing_structures.tsv";


    const string summary_file =
        output_root +
        "/structure_mapping_summary.txt";


    // ========================================================
    // STEP 1
    // 扫描所有PDB
    //
    // 建立：
    //
    // Viro3D ID -> PDB path
    // ========================================================

    cout
        << "========================================"
        << endl;

    cout
        << "Scanning Viro3D PDB structures..."
        << endl;

    cout
        << "========================================"
        << endl;


    unordered_map<string, string>
        viro3d_to_pdb;


    size_t total_pdb_files = 0;
    size_t valid_named_pdb = 0;
    size_t duplicate_pdb_id = 0;


    if (!fs::exists(pdb_root)) {

        cerr
            << "ERROR: PDB directory does not exist:\n"
            << pdb_root
            << endl;

        return 1;
    }


    for (
        const auto& entry :
        fs::recursive_directory_iterator(
            pdb_root
        )
    ) {

        if (!entry.is_regular_file()) {
            continue;
        }


        fs::path path =
            entry.path();


        string filename =
            path.filename().string();


        if (!endsWith(
                filename,
                ".pdb"
            )) {

            continue;
        }


        total_pdb_files++;


        string viro3d_id =
            pdbFilenameToViro3DID(
                filename
            );


        if (viro3d_id.empty()) {

            continue;
        }


        valid_named_pdb++;


        auto it =
            viro3d_to_pdb.find(
                viro3d_id
            );


        if (
            it !=
            viro3d_to_pdb.end()
        ) {

            duplicate_pdb_id++;

            cerr
                << "WARNING: duplicate Viro3D ID in PDB files: "
                << viro3d_id
                << endl;

        } else {

            viro3d_to_pdb[
                viro3d_id
            ] = path.string();
        }
    }


    cout
        << "Total PDB files:              "
        << total_pdb_files
        << endl;

    cout
        << "Valid *_relaxed.pdb files:    "
        << valid_named_pdb
        << endl;

    cout
        << "Unique structure IDs:         "
        << viro3d_to_pdb.size()
        << endl;

    cout
        << "Duplicate structure IDs:      "
        << duplicate_pdb_id
        << endl;


    // ========================================================
    // STEP 2
    // 读取 Baculoviridae CSV
    // ========================================================

    ifstream in(
        protein_csv
    );


    if (!in.is_open()) {

        cerr
            << "ERROR: cannot open protein CSV:\n"
            << protein_csv
            << endl;

        return 1;
    }


    string line;


    if (!getline(in, line)) {

        cerr
            << "ERROR: protein CSV is empty."
            << endl;

        return 1;
    }


    if (!line.empty() &&
        line.back() == '\r') {

        line.pop_back();
    }


    vector<string> headers =
        parseCSVLine(
            line
        );


    unordered_map<string, size_t>
        column_index;


    for (
        size_t i = 0;
        i < headers.size();
        ++i
    ) {

        column_index[
            trim(headers[i])
        ] = i;
    }


    // ========================================================
    // 检查关键列
    // ========================================================

    vector<string> required_columns = {

        "Viro3D ID",
        "ICTV Species",
        "Viro3D Name",
        "GenBank Protein ID",
        "UniProt ID",
        "Protein Length",
        "GenBank Genome Accession",
        "ColabFold pLDDT",
        "ColabFold pTM"
    };


    for (
        const string& col :
        required_columns
    ) {

        if (
            column_index.find(col)
            ==
            column_index.end()
        ) {

            cerr
                << "ERROR: missing column: "
                << col
                << endl;

            return 1;
        }
    }


    size_t idx_id =
        column_index["Viro3D ID"];

    size_t idx_species =
        column_index["ICTV Species"];

    size_t idx_name =
        column_index["Viro3D Name"];

    size_t idx_genbank_protein =
        column_index["GenBank Protein ID"];

    size_t idx_uniprot =
        column_index["UniProt ID"];

    size_t idx_length =
        column_index["Protein Length"];

    size_t idx_genome =
        column_index["GenBank Genome Accession"];

    size_t idx_plddt =
        column_index["ColabFold pLDDT"];

    size_t idx_ptm =
        column_index["ColabFold pTM"];


    // ========================================================
    // STEP 3
    // 打开输出文件
    // ========================================================

    ofstream mapping_out(
        mapping_file
    );


    ofstream missing_out(
        missing_file
    );


    if (
        !mapping_out.is_open() ||
        !missing_out.is_open()
    ) {

        cerr
            << "ERROR: cannot create output TSV files."
            << endl;

        return 1;
    }


    // ========================================================
    // mapping表
    // ========================================================

    mapping_out
        << "Viro3D_ID"
        << '\t'
        << "GenBank_Protein_ID"
        << '\t'
        << "UniProt_ID"
        << '\t'
        << "ICTV_Species"
        << '\t'
        << "Viro3D_Name"
        << '\t'
        << "Protein_Length"
        << '\t'
        << "Genome_Accession"
        << '\t'
        << "ColabFold_pLDDT"
        << '\t'
        << "ColabFold_pTM"
        << '\t'
        << "Structure_Status"
        << '\t'
        << "Original_PDB_Path"
        << '\t'
        << "Dataset_PDB_Path"
        << '\n';


    // ========================================================
    // missing表
    // ========================================================

    missing_out
        << "Viro3D_ID"
        << '\t'
        << "GenBank_Protein_ID"
        << '\t'
        << "ICTV_Species"
        << '\t'
        << "Viro3D_Name"
        << '\t'
        << "Protein_Length"
        << '\n';


    // ========================================================
    // 统计
    // ========================================================

    size_t total_proteins = 0;
    size_t matched = 0;
    size_t missing = 0;
    size_t symlink_success = 0;
    size_t symlink_failed = 0;


    unordered_set<string>
        seen_viro3d_ids;


    size_t duplicate_metadata_ids = 0;


    // ========================================================
    // STEP 4
    // 开始逐条匹配
    // ========================================================

    while (getline(in, line)) {

        if (!line.empty() &&
            line.back() == '\r') {

            line.pop_back();
        }


        if (line.empty()) {
            continue;
        }


        vector<string> fields =
            parseCSVLine(
                line
            );


        if (
            fields.size()
            <
            headers.size()
        ) {

            cerr
                << "WARNING: malformed CSV row skipped."
                << endl;

            continue;
        }


        total_proteins++;


        string viro3d_id =
            trim(
                fields[idx_id]
            );


        string species =
            trim(
                fields[idx_species]
            );


        string name =
            trim(
                fields[idx_name]
            );


        string genbank_id =
            trim(
                fields[idx_genbank_protein]
            );


        string uniprot =
            trim(
                fields[idx_uniprot]
            );


        string length =
            trim(
                fields[idx_length]
            );


        string genome =
            trim(
                fields[idx_genome]
            );


        string plddt =
            trim(
                fields[idx_plddt]
            );


        string ptm =
            trim(
                fields[idx_ptm]
            );


        // ====================================================
        // metadata ID重复检查
        // ====================================================

        if (
            !seen_viro3d_ids.insert(
                viro3d_id
            ).second
        ) {

            duplicate_metadata_ids++;
        }


        // ====================================================
        // 根据Viro3D ID查找PDB
        // ====================================================

        auto it =
            viro3d_to_pdb.find(
                viro3d_id
            );


        if (
            it !=
            viro3d_to_pdb.end()
        ) {

            // ================================================
            // MATCH
            // ================================================

            matched++;


            fs::path original_pdb =
                it->second;


            fs::path linked_pdb =
                fs::path(structure_dir)
                /
                original_pdb.filename();


            bool link_ok =
                createSymlinkSafe(
                    original_pdb,
                    linked_pdb
                );


            if (link_ok) {

                symlink_success++;

            } else {

                symlink_failed++;
            }


            mapping_out
                << viro3d_id
                << '\t'
                << genbank_id
                << '\t'
                << uniprot
                << '\t'
                << species
                << '\t'
                << name
                << '\t'
                << length
                << '\t'
                << genome
                << '\t'
                << plddt
                << '\t'
                << ptm
                << '\t'
                << "MATCHED"
                << '\t'
                << original_pdb.string()
                << '\t'
                << linked_pdb.string()
                << '\n';

        } else {

            // ================================================
            // MISSING
            // ================================================

            missing++;


            mapping_out
                << viro3d_id
                << '\t'
                << genbank_id
                << '\t'
                << uniprot
                << '\t'
                << species
                << '\t'
                << name
                << '\t'
                << length
                << '\t'
                << genome
                << '\t'
                << plddt
                << '\t'
                << ptm
                << '\t'
                << "MISSING_STRUCTURE"
                << '\t'
                << "NA"
                << '\t'
                << "NA"
                << '\n';


            missing_out
                << viro3d_id
                << '\t'
                << genbank_id
                << '\t'
                << species
                << '\t'
                << name
                << '\t'
                << length
                << '\n';
        }
    }


    in.close();

    mapping_out.close();

    missing_out.close();


    // ========================================================
    // STEP 5
    // summary
    // ========================================================

    double coverage = 0.0;


    if (total_proteins > 0) {

        coverage =
            100.0
            *
            static_cast<double>(matched)
            /
            static_cast<double>(
                total_proteins
            );
    }


    ofstream summary_out(
        summary_file
    );


    summary_out
        << "Baculoviridae structure mapping summary\n"
        << "========================================\n\n";


    summary_out
        << "Protein CSV:\n"
        << protein_csv
        << "\n\n";


    summary_out
        << "PDB source directory:\n"
        << pdb_root
        << "\n\n";


    summary_out
        << "Total Viro3D PDB files: "
        << total_pdb_files
        << "\n";


    summary_out
        << "Unique Viro3D structure IDs: "
        << viro3d_to_pdb.size()
        << "\n\n";


    summary_out
        << "Total Baculoviridae proteins: "
        << total_proteins
        << "\n";


    summary_out
        << "Matched structures: "
        << matched
        << "\n";


    summary_out
        << "Missing structures: "
        << missing
        << "\n";


    summary_out
        << fixed
        << setprecision(2);


    summary_out
        << "Structure coverage: "
        << coverage
        << "%\n\n";


    summary_out
        << "Symlinks created: "
        << symlink_success
        << "\n";


    summary_out
        << "Symlink failures: "
        << symlink_failed
        << "\n";


    summary_out
        << "Duplicate metadata Viro3D IDs: "
        << duplicate_metadata_ids
        << "\n";


    summary_out
        << "Duplicate PDB Viro3D IDs: "
        << duplicate_pdb_id
        << "\n";


    summary_out.close();


    // ========================================================
    // 打印结果
    // ========================================================

    cout
        << endl
        << "========================================"
        << endl;

    cout
        << "Baculoviridae structure mapping"
        << endl;

    cout
        << "========================================"
        << endl;


    cout
        << "Total proteins:           "
        << total_proteins
        << endl;

    cout
        << "Matched structures:       "
        << matched
        << endl;

    cout
        << "Missing structures:       "
        << missing
        << endl;


    cout
        << fixed
        << setprecision(2);


    cout
        << "Structure coverage:       "
        << coverage
        << "%"
        << endl;


    cout
        << "Symlinks created:         "
        << symlink_success
        << endl;

    cout
        << "Symlink failures:         "
        << symlink_failed
        << endl;

    cout
        << "Duplicate metadata IDs:   "
        << duplicate_metadata_ids
        << endl;

    cout
        << endl
        << "Structure directory:"
        << endl
        << structure_dir
        << endl;

    cout
        << endl
        << "Mapping:"
        << endl
        << mapping_file
        << endl;

    cout
        << endl
        << "Missing:"
        << endl
        << missing_file
        << endl;

    cout
        << endl
        << "Summary:"
        << endl
        << summary_file
        << endl;


    return 0;
}