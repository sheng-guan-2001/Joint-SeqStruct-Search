#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
namespace fs = std::filesystem;


// ============================================================
// 去掉字符串两端空白
// ============================================================
string trim(const string& s) {

    size_t start = 0;

    while (
        start < s.size() &&
        isspace(static_cast<unsigned char>(s[start]))
    ) {
        ++start;
    }

    size_t end = s.size();

    while (
        end > start &&
        isspace(static_cast<unsigned char>(s[end - 1]))
    ) {
        --end;
    }

    return s.substr(start, end - start);
}


// ============================================================
// 转为小写
// ============================================================
string toLower(string s) {

    transform(
        s.begin(),
        s.end(),
        s.begin(),
        [](unsigned char c) {
            return static_cast<char>(tolower(c));
        }
    );

    return s;
}


// ============================================================
// 合并连续空格
//
// "dna    polymerase"
// ->
// "dna polymerase"
// ============================================================
string collapseSpaces(const string& s) {

    string result;

    bool previous_space = false;

    for (unsigned char c : s) {

        if (isspace(c)) {

            if (!previous_space) {
                result += ' ';
                previous_space = true;
            }

        } else {

            result += static_cast<char>(c);
            previous_space = false;
        }
    }

    return trim(result);
}


// ============================================================
// 判断字符串是否以 prefix 开头
// ============================================================
bool startsWith(
    const string& text,
    const string& prefix
) {

    if (text.size() < prefix.size()) {
        return false;
    }

    return text.compare(
        0,
        prefix.size(),
        prefix
    ) == 0;
}


// ============================================================
// 去掉末尾无意义符号
//
// 例如：
// "p10;"
// ->
// "p10"
// ============================================================
string removeTrailingPunctuation(string s) {

    s = trim(s);

    while (!s.empty()) {

        char c = s.back();

        if (
            c == ';' ||
            c == ',' ||
            c == ':' ||
            c == '.' ||
            isspace(static_cast<unsigned char>(c))
        ) {

            s.pop_back();

        } else {

            break;
        }
    }

    return trim(s);
}


// ============================================================
// Viro3D Name 标准化
//
// 目前采用保守策略：
//
// 1. 去空格
// 2. 转小写
// 3. 去 "Product:" 前缀
// 4. 合并连续空格
// 5. 去末尾无意义标点
//
// 例如：
//
// Product: P10
// -> p10
//
// Product: ME53
// -> me53
//
// Product: DNA polymerase
// -> dna polymerase
//
// ============================================================
string normalizeViro3DName(
    const string& raw_name
) {

    string name =
        trim(raw_name);


    // 统一为小写
    name =
        toLower(name);


    // 合并连续空格
    name =
        collapseSpaces(name);


    // ========================================================
    // 去除 Product: 前缀
    // ========================================================

    const string product_prefix =
        "product:";


    if (
        startsWith(
            name,
            product_prefix
        )
    ) {

        name =
            name.substr(
                product_prefix.size()
            );

        name =
            trim(name);
    }


    // 再次规范空格
    name =
        collapseSpaces(name);


    // 去掉末尾标点
    name =
        removeTrailingPunctuation(
            name
        );


    // ========================================================
    // 如果整个字段为空
    // ========================================================

    if (name.empty()) {

        name =
            "unknown";
    }


    return name;
}


// ============================================================
// CSV parser
// 支持带逗号的双引号字段
// ============================================================
vector<string> parseCSVLine(
    const string& line
) {

    vector<string> fields;

    string current;

    bool in_quotes = false;


    for (
        size_t i = 0;
        i < line.size();
        ++i
    ) {

        char c =
            line[i];


        if (c == '"') {

            if (
                in_quotes &&
                i + 1 < line.size() &&
                line[i + 1] == '"'
            ) {

                current += '"';
                ++i;

            } else {

                in_quotes =
                    !in_quotes;
            }

        } else if (
            c == ',' &&
            !in_quotes
        ) {

            fields.push_back(
                current
            );

            current.clear();

        } else {

            current += c;
        }
    }


    fields.push_back(
        current
    );


    return fields;
}


// ============================================================
// CSV字段转义
// ============================================================
string escapeCSVField(
    const string& field
) {

    bool need_quotes =
        false;


    for (char c : field) {

        if (
            c == ',' ||
            c == '"' ||
            c == '\n' ||
            c == '\r'
        ) {

            need_quotes =
                true;

            break;
        }
    }


    if (!need_quotes) {

        return field;
    }


    string result =
        "\"";


    for (char c : field) {

        if (c == '"') {

            result +=
                "\"\"";

        } else {

            result += c;
        }
    }


    result +=
        "\"";


    return result;
}


// ============================================================
// 写CSV行
// ============================================================
void writeCSVRow(
    ofstream& out,
    const vector<string>& fields
) {

    for (
        size_t i = 0;
        i < fields.size();
        ++i
    ) {

        if (i > 0) {

            out << ",";
        }


        out
            << escapeCSVField(
                fields[i]
            );
    }


    out << "\n";
}


// ============================================================
// 主函数
// ============================================================
int main() {

    // ========================================================
    // 输入
    // ========================================================

    const string input_file =
        "Baculoviridae_proteins.csv";


    // ========================================================
    // 输出目录
    // ========================================================

    const string output_dir =
        "processed/Baculoviridae";


    fs::create_directories(
        output_dir
    );


    const string output_csv =
        output_dir +
        "/Baculoviridae_proteins_name_normalized.csv";


    const string mapping_file =
        output_dir +
        "/Viro3D_name_normalization_mapping.tsv";


    const string summary_file =
        output_dir +
        "/Viro3D_name_normalization_summary.txt";


    // ========================================================
    // 打开输入
    // ========================================================

    ifstream in(
        input_file
    );


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
            << "ERROR: empty CSV."
            << endl;

        return 1;
    }


    if (
        !line.empty() &&
        line.back() == '\r'
    ) {

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
    // 查找 Viro3D Name
    // ========================================================

    if (
        column_index.find(
            "Viro3D Name"
        )
        ==
        column_index.end()
    ) {

        cerr
            << "ERROR: Viro3D Name column not found."
            << endl;

        return 1;
    }


    size_t idx_name =
        column_index[
            "Viro3D Name"
        ];


    // ========================================================
    // 输出CSV
    // ========================================================

    ofstream csv_out(
        output_csv
    );


    if (!csv_out.is_open()) {

        cerr
            << "ERROR: cannot create output CSV."
            << endl;

        return 1;
    }


    // 原始列全部保留
    vector<string> new_headers =
        headers;


    new_headers.push_back(
        "Standardized Viro3D Name"
    );


    writeCSVRow(
        csv_out,
        new_headers
    );


    // ========================================================
    // 统计
    // ========================================================

    size_t total_rows =
        0;


    set<string>
        raw_names;


    set<string>
        normalized_names;


    // normalized name
    // ->
    // 原始names
    map<string, set<string>>
        normalized_to_raw;


    // raw name
    // ->
    // 出现次数
    map<string, size_t>
        raw_count;


    // normalized name
    // ->
    // 总protein数量
    map<string, size_t>
        normalized_count;


    // ========================================================
    // 逐行处理
    // ========================================================

    while (
        getline(
            in,
            line
        )
    ) {

        if (
            !line.empty() &&
            line.back() == '\r'
        ) {

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
                << "WARNING: malformed row skipped."
                << endl;

            continue;
        }


        total_rows++;


        string raw_name =
            trim(
                fields[
                    idx_name
                ]
            );


        string normalized_name =
            normalizeViro3DName(
                raw_name
            );


        raw_names.insert(
            raw_name
        );


        normalized_names.insert(
            normalized_name
        );


        normalized_to_raw[
            normalized_name
        ].insert(
            raw_name
        );


        raw_count[
            raw_name
        ]++;


        normalized_count[
            normalized_name
        ]++;


        // ====================================================
        // 新增标准化名称
        // ====================================================

        fields.push_back(
            normalized_name
        );


        writeCSVRow(
            csv_out,
            fields
        );
    }


    in.close();

    csv_out.close();


    // ========================================================
    // 输出 mapping
    // ========================================================

    ofstream mapping_out(
        mapping_file
    );


    mapping_out
        << "Standardized_Name"
        << '\t'
        << "Total_Protein_Count"
        << '\t'
        << "Raw_Name_Count"
        << '\t'
        << "Raw_Names"
        << '\n';


    for (
        const auto& pair :
        normalized_to_raw
    ) {

        const string&
            normalized =
            pair.first;


        const set<string>&
            raw_set =
            pair.second;


        mapping_out
            << normalized
            << '\t'

            << normalized_count[
                normalized
            ]
            << '\t'

            << raw_set.size()
            << '\t';


        bool first =
            true;


        for (
            const string& raw :
            raw_set
        ) {

            if (!first) {

                mapping_out
                    << " | ";
            }


            mapping_out
                << raw;


            first =
                false;
        }


        mapping_out
            << '\n';
    }


    mapping_out.close();


    // ========================================================
    // 统计有多少标准化名称合并了多个原始名称
    // ========================================================

    size_t merged_groups =
        0;


    size_t proteins_in_merged_groups =
        0;


    for (
        const auto& pair :
        normalized_to_raw
    ) {

        if (
            pair.second.size()
            > 1
        ) {

            merged_groups++;


            proteins_in_merged_groups +=
                normalized_count[
                    pair.first
                ];
        }
    }


    // ========================================================
    // summary
    // ========================================================

    ofstream summary_out(
        summary_file
    );


    summary_out
        << "Viro3D Name normalization summary\n";

    summary_out
        << "========================================\n\n";


    summary_out
        << "Total protein rows: "
        << total_rows
        << "\n";


    summary_out
        << "Unique raw Viro3D Names: "
        << raw_names.size()
        << "\n";


    summary_out
        << "Unique standardized Names: "
        << normalized_names.size()
        << "\n";


    summary_out
        << "Number of standardized groups containing >1 raw name: "
        << merged_groups
        << "\n";


    summary_out
        << "Proteins belonging to merged groups: "
        << proteins_in_merged_groups
        << "\n\n";


    summary_out
        << "Normalization rules:\n";

    summary_out
        << "  1. Convert to lowercase\n";

    summary_out
        << "  2. Remove leading 'Product:' prefix\n";

    summary_out
        << "  3. Trim leading/trailing spaces\n";

    summary_out
        << "  4. Collapse repeated whitespace\n";

    summary_out
        << "  5. Remove trailing punctuation\n";

    summary_out
        << "  6. No biological synonym merging is performed\n";


    summary_out.close();


    // ========================================================
    // 终端输出
    // ========================================================

    cout
        << "========================================"
        << endl;

    cout
        << "Viro3D Name normalization"
        << endl;

    cout
        << "========================================"
        << endl;


    cout
        << "Total proteins:                  "
        << total_rows
        << endl;


    cout
        << "Unique raw names:                "
        << raw_names.size()
        << endl;


    cout
        << "Unique standardized names:       "
        << normalized_names.size()
        << endl;


    cout
        << "Merged standardized groups:      "
        << merged_groups
        << endl;


    cout
        << "Proteins in merged groups:       "
        << proteins_in_merged_groups
        << endl;


    cout
        << endl
        << "Output normalized CSV:"
        << endl
        << output_csv
        << endl;


    cout
        << endl
        << "Name mapping:"
        << endl
        << mapping_file
        << endl;


    cout
        << endl
        << "Summary:"
        << endl
        << summary_file
        << endl;


    return 0;
}