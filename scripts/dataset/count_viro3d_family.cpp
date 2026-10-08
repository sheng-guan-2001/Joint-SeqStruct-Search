#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

//找出family的代码


// ============================================================
// 去除字符串两端空白
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
// 支持字段内部带逗号、双引号
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
// Family统计结构
// ============================================================
struct FamilyStats {

    // 不同 species
    set<string> species;

    // 不同 genome accession
    set<string> genomes;

    // 蛋白条数
    size_t protein_count = 0;
};


// ============================================================
// 输出结构
// ============================================================
struct OutputRow {

    string family;

    size_t species_count;

    size_t genome_count;

    size_t protein_count;

    double proteins_per_species;

    double proteins_per_genome;
};


// ============================================================
// 主函数
// ============================================================
int main() {

    // ========================================================
    // 输入文件
    // ========================================================

    const string virus_file =
        "Viro3D_viruses_list_with_coverage.csv";


    const string protein_file =
        "Viro3D_proteins_filtered.csv";


    // ========================================================
    // 输出文件
    // ========================================================

    const string output_dir =
        "protein_count_by_ICTV_family";


    fs::create_directories(output_dir);


    const string output_file =
        output_dir +
        "/protein_count_by_ICTV_family.tsv";


    // ========================================================
    // 第一部分：
    // 读取 virus metadata
    //
    // 建立：
    //
    // ICTV Sort -> ICTV Family
    //
    // 例如：
    // 33 -> Alloherpesviridae
    // ========================================================

    ifstream virus_in(virus_file);

    if (!virus_in.is_open()) {

        cerr
            << "ERROR: cannot open virus file:\n"
            << virus_file
            << endl;

        return 1;
    }


    string header_line;

    if (!getline(virus_in, header_line)) {

        cerr
            << "ERROR: virus CSV is empty."
            << endl;

        return 1;
    }


    if (!header_line.empty() &&
        header_line.back() == '\r') {

        header_line.pop_back();
    }


    vector<string> virus_headers =
        parseCSVLine(header_line);


    unordered_map<string, size_t>
        virus_col;


    for (size_t i = 0;
         i < virus_headers.size();
         ++i) {

        virus_col[
            trim(virus_headers[i])
        ] = i;
    }


    // ========================================================
    // 检查列
    // ========================================================

    vector<string> required_virus_columns = {

        "ICTV Sort",
        "ICTV Family",
        "ICTV Species"
    };


    for (const auto& col :
         required_virus_columns) {

        if (
            virus_col.find(col)
            ==
            virus_col.end()
        ) {

            cerr
                << "ERROR: missing column in virus CSV: "
                << col
                << endl;

            return 1;
        }
    }


    size_t idx_v_sort =
        virus_col["ICTV Sort"];

    size_t idx_v_family =
        virus_col["ICTV Family"];

    size_t idx_v_species =
        virus_col["ICTV Species"];


    // ========================================================
    // ICTV Sort -> Family
    // ========================================================

    unordered_map<string, string>
        sort_to_family;


    // 同时保存 species -> family
    // 作为备用检查
    unordered_map<string, string>
        species_to_family;


    string line;

    size_t virus_rows = 0;


    while (getline(virus_in, line)) {

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
            virus_headers.size()
        ) {

            cerr
                << "WARNING: malformed virus row skipped."
                << endl;

            continue;
        }


        string sort_id =
            trim(fields[idx_v_sort]);

        string family =
            trim(fields[idx_v_family]);

        string species =
            trim(fields[idx_v_species]);


        if (family.empty()) {
            family = "UNKNOWN_FAMILY";
        }


        if (!sort_id.empty()) {

            sort_to_family[
                sort_id
            ] = family;
        }


        if (!species.empty()) {

            species_to_family[
                species
            ] = family;
        }


        virus_rows++;
    }


    virus_in.close();


    cout
        << "Virus metadata rows loaded: "
        << virus_rows
        << endl;

    cout
        << "ICTV Sort -> Family mappings: "
        << sort_to_family.size()
        << endl;


    // ========================================================
    // 第二部分：
    // 读取已经过滤后的 protein CSV
    // ========================================================

    ifstream protein_in(
        protein_file
    );


    if (!protein_in.is_open()) {

        cerr
            << "ERROR: cannot open protein file:\n"
            << protein_file
            << endl;

        return 1;
    }


    if (!getline(
            protein_in,
            header_line
        )) {

        cerr
            << "ERROR: protein CSV is empty."
            << endl;

        return 1;
    }


    if (!header_line.empty() &&
        header_line.back() == '\r') {

        header_line.pop_back();
    }


    vector<string> protein_headers =
        parseCSVLine(header_line);


    unordered_map<string, size_t>
        protein_col;


    for (size_t i = 0;
         i < protein_headers.size();
         ++i) {

        protein_col[
            trim(
                protein_headers[i]
            )
        ] = i;
    }


    // ========================================================
    // 检查 protein 文件关键字段
    // ========================================================

    vector<string> required_protein_columns = {

        "ICTV Sort",
        "ICTV Species",
        "GenBank Genome Accession"
    };


    for (
        const auto& col :
        required_protein_columns
    ) {

        if (
            protein_col.find(col)
            ==
            protein_col.end()
        ) {

            cerr
                << "ERROR: missing column in protein CSV: "
                << col
                << endl;

            return 1;
        }
    }


    size_t idx_p_sort =
        protein_col[
            "ICTV Sort"
        ];

    size_t idx_p_species =
        protein_col[
            "ICTV Species"
        ];

    size_t idx_p_genome =
        protein_col[
            "GenBank Genome Accession"
        ];


    // ========================================================
    // Family -> 统计结果
    // ========================================================

    unordered_map<string, FamilyStats>
        family_stats;


    size_t total_proteins = 0;

    size_t matched_family = 0;

    size_t unmatched_family = 0;


    // ========================================================
    // 遍历每条 protein
    // ========================================================

    while (getline(
        protein_in,
        line
    )) {

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
            protein_headers.size()
        ) {

            cerr
                << "WARNING: malformed protein row skipped."
                << endl;

            continue;
        }


        total_proteins++;


        string sort_id =
            trim(
                fields[idx_p_sort]
            );

        string species =
            trim(
                fields[idx_p_species]
            );

        string genome =
            trim(
                fields[idx_p_genome]
            );


        // ====================================================
        // 找 family
        //
        // 优先使用 ICTV Sort
        // 如果 ICTV Sort 找不到，
        // 再用 species 备用匹配
        // ====================================================

        string family;


        auto it_sort =
            sort_to_family.find(
                sort_id
            );


        if (
            it_sort
            !=
            sort_to_family.end()
        ) {

            family =
                it_sort->second;

            matched_family++;

        } else {

            auto it_species =
                species_to_family.find(
                    species
                );


            if (
                it_species
                !=
                species_to_family.end()
            ) {

                family =
                    it_species->second;

                matched_family++;

            } else {

                family =
                    "UNKNOWN_FAMILY";

                unmatched_family++;
            }
        }


        // ====================================================
        // 加入该 family 统计
        // ====================================================

        FamilyStats&
            stats =
            family_stats[family];


        // 每一行就是一个 protein
        stats.protein_count++;


        // set 自动去重
        if (!species.empty()) {

            stats.species.insert(
                species
            );
        }


        if (!genome.empty()) {

            stats.genomes.insert(
                genome
            );
        }
    }


    protein_in.close();


    // ========================================================
    // 转为 vector 方便排序
    // ========================================================

    vector<OutputRow>
        rows;


    for (
        const auto& pair :
        family_stats
    ) {

        const string& family =
            pair.first;

        const FamilyStats& stats =
            pair.second;


        OutputRow row;


        row.family =
            family;


        row.species_count =
            stats.species.size();


        row.genome_count =
            stats.genomes.size();


        row.protein_count =
            stats.protein_count;


        if (
            row.species_count > 0
        ) {

            row.proteins_per_species =
                static_cast<double>(
                    row.protein_count
                )
                /
                static_cast<double>(
                    row.species_count
                );

        } else {

            row.proteins_per_species =
                0.0;
        }


        if (
            row.genome_count > 0
        ) {

            row.proteins_per_genome =
                static_cast<double>(
                    row.protein_count
                )
                /
                static_cast<double>(
                    row.genome_count
                );

        } else {

            row.proteins_per_genome =
                0.0;
        }


        rows.push_back(row);
    }


    // ========================================================
    // 按 Protein Count 从高到低
    // ========================================================

    sort(
        rows.begin(),
        rows.end(),

        [](
            const OutputRow& a,
            const OutputRow& b
        ) {

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
                a.family
                <
                b.family;
        }
    );


    // ========================================================
    // 写输出
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
        << "ICTV_Family"
        << '\t'
        << "Species_Count"
        << '\t'
        << "Genome_Count"
        << '\t'
        << "Protein_Count"
        << '\t'
        << "Mean_Proteins_Per_Species"
        << '\t'
        << "Mean_Proteins_Per_Genome"
        << '\n';


    out
        << fixed
        << setprecision(2);


    for (
        const auto& row :
        rows
    ) {

        out
            << row.family
            << '\t'

            << row.species_count
            << '\t'

            << row.genome_count
            << '\t'

            << row.protein_count
            << '\t'

            << row.proteins_per_species
            << '\t'

            << row.proteins_per_genome
            << '\n';
    }


    out.close();


    // ========================================================
    // 终端输出 summary
    // ========================================================

    cout
        << endl
        << "========================================"
        << endl;

    cout
        << "ICTV Family summary"
        << endl;

    cout
        << "========================================"
        << endl;


    cout
        << "Filtered proteins:          "
        << total_proteins
        << endl;

    cout
        << "Family-matched proteins:    "
        << matched_family
        << endl;

    cout
        << "Unmatched proteins:         "
        << unmatched_family
        << endl;

    cout
        << "Unique ICTV families:       "
        << family_stats.size()
        << endl;


    cout
        << endl
        << "Top 30 families by protein count:"
        << endl;

    cout
        << "----------------------------------------"
        << endl;


    cout
        << left
        << setw(35)
        << "ICTV Family"

        << right
        << setw(12)
        << "Species"

        << setw(12)
        << "Genomes"

        << setw(12)
        << "Proteins"

        << endl;


    cout
        << string(
            71,
            '-'
        )
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
            << left
            << setw(35)
            << rows[i].family

            << right
            << setw(12)
            << rows[i].species_count

            << setw(12)
            << rows[i].genome_count

            << setw(12)
            << rows[i].protein_count

            << endl;
    }


    cout
        << endl
        << "Output file:"
        << endl
        << output_file
        << endl;


    return 0;
}