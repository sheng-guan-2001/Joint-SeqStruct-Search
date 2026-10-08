#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
namespace fs = std::filesystem;


// ============================================================
// 去除字符串两端空格
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
// 支持：
// abc,def,"hello,world",123
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
// 安全转换 double
// ============================================================
bool parseDouble(
    const string& text,
    double& value
) {

    string x = trim(text);

    if (x.empty()) {
        return false;
    }

    try {

        size_t pos = 0;

        value = stod(x, &pos);

        return pos == x.size();

    } catch (...) {

        return false;
    }
}


// ============================================================
// 安全转换整数
// ============================================================
bool parseLongLong(
    const string& text,
    long long& value
) {

    string x = trim(text);

    if (x.empty()) {
        return false;
    }

    try {

        size_t pos = 0;

        value = stoll(x, &pos);

        return pos == x.size();

    } catch (...) {

        return false;
    }
}


// ============================================================
// 每一个 Viro3D Name 的统计结构
// ============================================================
struct NameStats {

    size_t protein_count = 0;

    set<string> species;
    set<string> genomes;

    double sum_plddt = 0.0;
    double sum_ptm = 0.0;
    double sum_length = 0.0;

    size_t valid_plddt = 0;
    size_t valid_ptm = 0;
    size_t valid_length = 0;
};


// ============================================================
// 输出结构
// ============================================================
struct OutputRow {

    string name;

    size_t protein_count = 0;
    size_t species_count = 0;
    size_t genome_count = 0;

    double mean_plddt = 0.0;
    double mean_ptm = 0.0;
    double mean_length = 0.0;
};


// ============================================================
// 主函数
// ============================================================
int main() {

    // ========================================================
    // 输入：
    // 已经筛选好的 Baculoviridae 数据
    // ========================================================

    const string input_file =
        "Baculoviridae_proteins_name_normalized.csv";


    // ========================================================
    // 输出目录
    // ========================================================

    const string output_dir =
        "processed/Baculoviridae";


    fs::create_directories(output_dir);


    const string output_file =
        output_dir +
        "/Viro3D_name_statistics.tsv";


    const string summary_file =
        output_dir +
        "/Viro3D_name_statistics_summary.txt";


    // ========================================================
    // 打开输入文件
    // ========================================================

    ifstream in(input_file);

    if (!in.is_open()) {

        cerr
            << "ERROR: cannot open input file:\n"
            << input_file
            << endl;

        return 1;
    }


    // ========================================================
    // 读取表头
    // ========================================================

    string line;

    if (!getline(in, line)) {

        cerr
            << "ERROR: input CSV is empty."
            << endl;

        return 1;
    }


    if (!line.empty() &&
        line.back() == '\r') {

        line.pop_back();
    }


    vector<string> headers =
        parseCSVLine(line);


    unordered_map<string, size_t>
        column_index;


    for (size_t i = 0;
         i < headers.size();
         ++i) {

        column_index[
            trim(headers[i])
        ] = i;
    }


    // ========================================================
    // 检查关键列
    // ========================================================

    vector<string> required_columns = {

        "Standardized Viro3D Name",
        "ICTV Species",
        "GenBank Genome Accession",
        "Protein Length",
        "ColabFold pLDDT",
        "ColabFold pTM"
    };


    for (const string& col :
         required_columns) {

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


    size_t idx_name =
        column_index["Standardized Viro3D Name"];

    size_t idx_species =
        column_index["ICTV Species"];

    size_t idx_genome =
        column_index["GenBank Genome Accession"];

    size_t idx_length =
        column_index["Protein Length"];

    size_t idx_plddt =
        column_index["ColabFold pLDDT"];

    size_t idx_ptm =
        column_index["ColabFold pTM"];


    // ========================================================
    // Viro3D Name -> stats
    // ========================================================

    unordered_map<string, NameStats>
        stats_map;


    size_t total_rows = 0;


    // ========================================================
    // 遍历所有蛋白
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
            parseCSVLine(line);


        if (
            fields.size()
            <
            headers.size()
        ) {

            cerr
                << "WARNING: malformed row skipped."
                << endl;

            continue;
        }


        total_rows++;


        string name =
            trim(fields[idx_name]);


        string species =
            trim(fields[idx_species]);


        string genome =
            trim(fields[idx_genome]);


        if (name.empty()) {

            name =
                "UNKNOWN_NAME";
        }


        NameStats& stats =
            stats_map[name];


        // ====================================================
        // Protein Count
        // ====================================================

        stats.protein_count++;


        // ====================================================
        // Species Count
        // ====================================================

        if (!species.empty()) {

            stats.species.insert(
                species
            );
        }


        // ====================================================
        // Genome Count
        // ====================================================

        if (!genome.empty()) {

            stats.genomes.insert(
                genome
            );
        }


        // ====================================================
        // Length
        // ====================================================

        long long length = 0;

        if (
            parseLongLong(
                fields[idx_length],
                length
            )
        ) {

            stats.sum_length +=
                static_cast<double>(
                    length
                );

            stats.valid_length++;
        }


        // ====================================================
        // pLDDT
        // ====================================================

        double plddt = 0.0;

        if (
            parseDouble(
                fields[idx_plddt],
                plddt
            )
        ) {

            stats.sum_plddt += plddt;

            stats.valid_plddt++;
        }


        // ====================================================
        // pTM
        // ====================================================

        double ptm = 0.0;

        if (
            parseDouble(
                fields[idx_ptm],
                ptm
            )
        ) {

            stats.sum_ptm += ptm;

            stats.valid_ptm++;
        }
    }


    in.close();


    // ========================================================
    // 转换为 vector
    // ========================================================

    vector<OutputRow>
        rows;


    for (
        const auto& pair :
        stats_map
    ) {

        const string& name =
            pair.first;

        const NameStats& stats =
            pair.second;


        OutputRow row;


        row.name =
            name;


        row.protein_count =
            stats.protein_count;


        row.species_count =
            stats.species.size();


        row.genome_count =
            stats.genomes.size();


        if (stats.valid_plddt > 0) {

            row.mean_plddt =
                stats.sum_plddt
                /
                static_cast<double>(
                    stats.valid_plddt
                );
        }


        if (stats.valid_ptm > 0) {

            row.mean_ptm =
                stats.sum_ptm
                /
                static_cast<double>(
                    stats.valid_ptm
                );
        }


        if (stats.valid_length > 0) {

            row.mean_length =
                stats.sum_length
                /
                static_cast<double>(
                    stats.valid_length
                );
        }


        rows.push_back(
            row
        );
    }


    // ========================================================
    // 排序
    //
    // 第一优先：
    // Species Count 高
    //
    // 第二优先：
    // Protein Count 高
    //
    // 这样特别适合后面挑 query
    // ========================================================

    sort(
        rows.begin(),
        rows.end(),

        [](
            const OutputRow& a,
            const OutputRow& b
        ) {

            if (
                a.species_count
                !=
                b.species_count
            ) {

                return
                    a.species_count
                    >
                    b.species_count;
            }


            if (
                a.protein_count
                !=
                b.protein_count
            ) {

                return
                    a.protein_count
                    >
                    b.protein_count;
            }


            return
                a.name
                <
                b.name;
        }
    );


    // ========================================================
    // 写 TSV
    // ========================================================

    ofstream out(
        output_file
    );


    if (!out.is_open()) {

        cerr
            << "ERROR: cannot create output file."
            << endl;

        return 1;
    }


    out
        << "Viro3D_Name"
        << '\t'
        << "Protein_Count"
        << '\t'
        << "Species_Count"
        << '\t'
        << "Genome_Count"
        << '\t'
        << "Mean_ColabFold_pLDDT"
        << '\t'
        << "Mean_ColabFold_pTM"
        << '\t'
        << "Mean_Protein_Length"
        << '\n';


    out
        << fixed
        << setprecision(2);


    for (
        const OutputRow& row :
        rows
    ) {

        out
            << row.name
            << '\t'

            << row.protein_count
            << '\t'

            << row.species_count
            << '\t'

            << row.genome_count
            << '\t'

            << row.mean_plddt
            << '\t'

            << row.mean_ptm
            << '\t'

            << row.mean_length
            << '\n';
    }


    out.close();


    // ========================================================
    // 简单summary
    // ========================================================

    size_t names_species_5 = 0;
    size_t names_species_10 = 0;
    size_t names_species_20 = 0;
    size_t names_species_50 = 0;


    for (
        const OutputRow& row :
        rows
    ) {

        if (row.species_count >= 5) {
            names_species_5++;
        }

        if (row.species_count >= 10) {
            names_species_10++;
        }

        if (row.species_count >= 20) {
            names_species_20++;
        }

        if (row.species_count >= 50) {
            names_species_50++;
        }
    }


    ofstream summary(
        summary_file
    );


    summary
        << "Baculoviridae Viro3D Name statistics\n";

    summary
        << "========================================\n\n";


    summary
        << "Total proteins: "
        << total_rows
        << "\n";


    summary
        << "Unique Viro3D Names: "
        << rows.size()
        << "\n\n";


    summary
        << "Names present in >= 5 species: "
        << names_species_5
        << "\n";


    summary
        << "Names present in >= 10 species: "
        << names_species_10
        << "\n";


    summary
        << "Names present in >= 20 species: "
        << names_species_20
        << "\n";


    summary
        << "Names present in >= 50 species: "
        << names_species_50
        << "\n";


    summary.close();


    // ========================================================
    // 终端输出
    // ========================================================

    cout
        << "========================================"
        << endl;

    cout
        << "Baculoviridae Viro3D Name statistics"
        << endl;

    cout
        << "========================================"
        << endl;


    cout
        << "Total proteins:        "
        << total_rows
        << endl;

    cout
        << "Unique Viro3D Names:   "
        << rows.size()
        << endl;

    cout
        << "Names in >=5 species:  "
        << names_species_5
        << endl;

    cout
        << "Names in >=10 species: "
        << names_species_10
        << endl;

    cout
        << "Names in >=20 species: "
        << names_species_20
        << endl;

    cout
        << "Names in >=50 species: "
        << names_species_50
        << endl;


    cout
        << endl
        << "Top 30 by Species Count:"
        << endl;

    cout
        << "----------------------------------------"
        << endl;


    size_t top_n =
        min<size_t>(
            30,
            rows.size()
        );


    for (
        size_t i = 0;
        i < top_n;
        ++i
    ) {

        cout
            << setw(3)
            << i + 1
            << "  "

            << setw(4)
            << rows[i].species_count
            << " species  "

            << setw(5)
            << rows[i].protein_count
            << " proteins  "

            << rows[i].name
            << endl;
    }


    cout
        << endl
        << "Output:"
        << endl
        << output_file
        << endl;


    cout
        << endl
        << "Summary:"
        << endl
        << summary_file
        << endl;


    return 0;
}