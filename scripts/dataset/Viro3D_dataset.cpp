#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;


//筛选protein的文件

// ============================================================
// 去掉字符串两端空白
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
// CSV解析器
//
// 支持：
// abc,def,"hello,world",123
//
// 也支持CSV内部双引号：
// "abc ""test"" xyz"
// ============================================================
vector<string> parseCSVLine(const string& line) {

    vector<string> fields;

    string current;

    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {

        char c = line[i];

        if (c == '"') {

            // CSV中两个连续引号表示真正的 "
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
// 将一行字段重新写成CSV
//
// 如果字段中含逗号、引号、换行，自动加引号
// ============================================================
string escapeCSVField(const string& field) {

    bool need_quotes = false;

    for (char c : field) {
        if (c == ',' || c == '"' ||
            c == '\n' || c == '\r') {
            need_quotes = true;
            break;
        }
    }

    if (!need_quotes) {
        return field;
    }

    string result = "\"";

    for (char c : field) {

        if (c == '"') {
            result += "\"\"";
        } else {
            result += c;
        }
    }

    result += "\"";

    return result;
}


// ============================================================
// 写CSV行
// ============================================================
void writeCSVRow(
    ofstream& out,
    const vector<string>& fields
) {

    for (size_t i = 0; i < fields.size(); ++i) {

        if (i > 0) {
            out << ",";
        }

        out << escapeCSVField(fields[i]);
    }

    out << "\n";
}


// ============================================================
// 安全转换 double
// ============================================================
bool parseDouble(
    const string& s,
    double& value
) {

    string x = trim(s);

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
bool parseLong(
    const string& s,
    long long& value
) {

    string x = trim(s);

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
// 主程序
// ============================================================
int main() {

    // ========================================================
    // 输入文件
    // ========================================================

    const string input_file =
        "Viro3D_proteins_list.csv";


    // ========================================================
    // 输出目录
    // ========================================================

    const string output_dir =
        "processed/filtered_proteins";


    fs::create_directories(output_dir);


    const string output_csv =
        output_dir +
        "/Viro3D_proteins_filtered.csv";


    const string summary_file =
        output_dir +
        "/filter_summary.txt";


    // ========================================================
    // 筛选阈值
    // ========================================================

    const string REQUIRED_CATEGORY =
        "protein";

    const double MIN_PLDDT =
        80.0;

    const long long MIN_LENGTH =
        50;

    const long long MAX_LENGTH =
        2000;


    // ========================================================
    // 打开文件
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

    string header_line;

    if (!getline(in, header_line)) {

        cerr
            << "ERROR: empty CSV file."
            << endl;

        return 1;
    }


    if (!header_line.empty() &&
        header_line.back() == '\r') {

        header_line.pop_back();
    }


    vector<string> headers =
        parseCSVLine(header_line);


    cout
        << "========================================"
        << endl;

    cout
        << "Viro3D protein filtering"
        << endl;

    cout
        << "========================================"
        << endl;


    cout
        << "Columns detected: "
        << headers.size()
        << endl;


    // ========================================================
    // 建立字段名 -> 列号
    // ========================================================

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
    // 检查必须存在的列
    // ========================================================

    vector<string> required_columns = {

        "Peptide Category",
        "Protein Length",
        "ColabFold pLDDT"
    };


    for (const string& col :
         required_columns) {

        if (
            column_index.find(col)
            ==
            column_index.end()
        ) {

            cerr
                << "ERROR: required column not found: "
                << col
                << endl;

            return 1;
        }
    }


    size_t idx_category =
        column_index["Peptide Category"];

    size_t idx_length =
        column_index["Protein Length"];

    size_t idx_plddt =
        column_index["ColabFold pLDDT"];


    // 可选列
    bool has_uniprot =
        column_index.find("UniProt ID")
        != column_index.end();

    bool has_ptm =
        column_index.find("ColabFold pTM")
        != column_index.end();

    bool has_sequence =
        column_index.find("Protein Sequence")
        != column_index.end();

    bool has_viro3d_id =
        column_index.find("Viro3D ID")
        != column_index.end();


    size_t idx_uniprot = 0;
    size_t idx_ptm = 0;
    size_t idx_sequence = 0;
    size_t idx_viro3d_id = 0;


    if (has_uniprot) {
        idx_uniprot =
            column_index["UniProt ID"];
    }

    if (has_ptm) {
        idx_ptm =
            column_index["ColabFold pTM"];
    }

    if (has_sequence) {
        idx_sequence =
            column_index["Protein Sequence"];
    }

    if (has_viro3d_id) {
        idx_viro3d_id =
            column_index["Viro3D ID"];
    }


    // ========================================================
    // 打开输出CSV
    // ========================================================

    ofstream out(output_csv);

    if (!out.is_open()) {

        cerr
            << "ERROR: cannot create output file:\n"
            << output_csv
            << endl;

        return 1;
    }


    // 写原始表头
    writeCSVRow(
        out,
        headers
    );


    // ========================================================
    // 统计变量
    // ========================================================

    size_t total_rows = 0;

    size_t protein_rows = 0;

    size_t protein_with_valid_length = 0;

    size_t protein_length_pass = 0;

    size_t protein_valid_plddt = 0;

    size_t final_pass = 0;


    size_t missing_length = 0;
    size_t missing_plddt = 0;

    size_t missing_uniprot_after_filter = 0;
    size_t missing_ptm_after_filter = 0;
    size_t missing_sequence_after_filter = 0;


    double sum_length = 0.0;
    double sum_plddt = 0.0;
    double sum_ptm = 0.0;

    size_t valid_ptm_after_filter = 0;


    long long min_selected_length =
        numeric_limits<long long>::max();

    long long max_selected_length =
        numeric_limits<long long>::min();

    double min_selected_plddt =
        numeric_limits<double>::max();

    double max_selected_plddt =
        numeric_limits<double>::lowest();


    // ========================================================
    // 逐行读取CSV
    // ========================================================

    string line;

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


        // 如果这一行字段数少于表头，
        // 说明可能存在格式异常
        if (fields.size() < headers.size()) {

            cerr
                << "WARNING: malformed row, "
                << "expected "
                << headers.size()
                << " columns but got "
                << fields.size()
                << endl;

            continue;
        }


        total_rows++;


        // ====================================================
        // Step 1:
        // Peptide Category == protein
        // ====================================================

        string category =
            trim(fields[idx_category]);


        if (category != REQUIRED_CATEGORY) {
            continue;
        }


        protein_rows++;


        // ====================================================
        // Step 2:
        // Protein Length有效
        // ====================================================

        long long protein_length = 0;


        if (
            !parseLong(
                fields[idx_length],
                protein_length
            )
        ) {

            missing_length++;
            continue;
        }


        protein_with_valid_length++;


        // ====================================================
        // Step 3:
        // 50 <= Protein Length <= 2000
        // ====================================================

        if (
            protein_length < MIN_LENGTH ||
            protein_length > MAX_LENGTH
        ) {

            continue;
        }


        protein_length_pass++;


        // ====================================================
        // Step 4:
        // pLDDT有效
        // ====================================================

        double plddt = 0.0;


        if (
            !parseDouble(
                fields[idx_plddt],
                plddt
            )
        ) {

            missing_plddt++;
            continue;
        }


        protein_valid_plddt++;


        // ====================================================
        // Step 5:
        // ColabFold pLDDT >= 80
        // ====================================================

        if (plddt < MIN_PLDDT) {
            continue;
        }


        // ====================================================
        // 最终通过
        // ====================================================

        final_pass++;


        // 写入完整CSV行
        writeCSVRow(
            out,
            fields
        );


        // ====================================================
        // 统计
        // ====================================================

        sum_length +=
            static_cast<double>(
                protein_length
            );

        sum_plddt += plddt;


        min_selected_length =
            min(
                min_selected_length,
                protein_length
            );


        max_selected_length =
            max(
                max_selected_length,
                protein_length
            );


        min_selected_plddt =
            min(
                min_selected_plddt,
                plddt
            );


        max_selected_plddt =
            max(
                max_selected_plddt,
                plddt
            );


        // ----------------------------------------------------
        // UniProt缺失统计
        // ----------------------------------------------------

        if (has_uniprot) {

            string uniprot =
                trim(
                    fields[
                        idx_uniprot
                    ]
                );


            if (
                uniprot.empty() ||
                uniprot == "NA" ||
                uniprot == "N/A"
            ) {

                missing_uniprot_after_filter++;
            }
        }


        // ----------------------------------------------------
        // pTM统计
        // ----------------------------------------------------

        if (has_ptm) {

            double ptm = 0.0;


            if (
                parseDouble(
                    fields[idx_ptm],
                    ptm
                )
            ) {

                sum_ptm += ptm;

                valid_ptm_after_filter++;

            } else {

                missing_ptm_after_filter++;
            }
        }


        // ----------------------------------------------------
        // sequence缺失统计
        // ----------------------------------------------------

        if (has_sequence) {

            string seq =
                trim(
                    fields[
                        idx_sequence
                    ]
                );


            if (seq.empty()) {

                missing_sequence_after_filter++;
            }
        }
    }


    in.close();
    out.close();


    // ========================================================
    // 计算均值
    // ========================================================

    double mean_length =
        final_pass > 0
        ? sum_length /
          static_cast<double>(
              final_pass
          )
        : 0.0;


    double mean_plddt =
        final_pass > 0
        ? sum_plddt /
          static_cast<double>(
              final_pass
          )
        : 0.0;


    double mean_ptm =
        valid_ptm_after_filter > 0
        ? sum_ptm /
          static_cast<double>(
              valid_ptm_after_filter
          )
        : 0.0;


    // ========================================================
    // 输出summary文件
    // ========================================================

    ofstream summary(
        summary_file
    );


    summary
        << "Viro3D protein filtering summary\n";

    summary
        << "========================================\n\n";


    summary
        << "Input file:\n"
        << input_file
        << "\n\n";


    summary
        << "Filtering criteria:\n";

    summary
        << "  Peptide Category == protein\n";

    summary
        << "  ColabFold pLDDT >= "
        << MIN_PLDDT
        << "\n";

    summary
        << "  Protein Length >= "
        << MIN_LENGTH
        << "\n";

    summary
        << "  Protein Length <= "
        << MAX_LENGTH
        << "\n\n";


    summary
        << "Counts:\n";

    summary
        << "  Total input rows: "
        << total_rows
        << "\n";

    summary
        << "  Category == protein: "
        << protein_rows
        << "\n";

    summary
        << "  Protein with valid length: "
        << protein_with_valid_length
        << "\n";

    summary
        << "  Protein passing length filter: "
        << protein_length_pass
        << "\n";

    summary
        << "  Length-filtered proteins with valid pLDDT: "
        << protein_valid_plddt
        << "\n";

    summary
        << "  Final selected proteins: "
        << final_pass
        << "\n\n";


    summary
        << "Missing values:\n";

    summary
        << "  Missing/invalid Protein Length: "
        << missing_length
        << "\n";

    summary
        << "  Missing/invalid ColabFold pLDDT: "
        << missing_plddt
        << "\n";

    if (has_uniprot) {

        summary
            << "  Missing UniProt ID after filtering: "
            << missing_uniprot_after_filter
            << "\n";
    }

    if (has_ptm) {

        summary
            << "  Missing ColabFold pTM after filtering: "
            << missing_ptm_after_filter
            << "\n";
    }

    if (has_sequence) {

        summary
            << "  Missing Protein Sequence after filtering: "
            << missing_sequence_after_filter
            << "\n";
    }


    if (final_pass > 0) {

        summary
            << "\nSelected dataset statistics:\n";

        summary
            << fixed
            << setprecision(2);

        summary
            << "  Mean protein length: "
            << mean_length
            << " aa\n";

        summary
            << "  Minimum protein length: "
            << min_selected_length
            << " aa\n";

        summary
            << "  Maximum protein length: "
            << max_selected_length
            << " aa\n";

        summary
            << "  Mean ColabFold pLDDT: "
            << mean_plddt
            << "\n";

        summary
            << "  Minimum ColabFold pLDDT: "
            << min_selected_plddt
            << "\n";

        summary
            << "  Maximum ColabFold pLDDT: "
            << max_selected_plddt
            << "\n";

        if (
            has_ptm &&
            valid_ptm_after_filter > 0
        ) {

            summary
                << "  Mean ColabFold pTM: "
                << mean_ptm
                << "\n";
        }
    }


    summary.close();


    // ========================================================
    // 输出到终端
    // ========================================================

    cout
        << endl
        << "========================================"
        << endl;

    cout
        << "Filtering summary"
        << endl;

    cout
        << "========================================"
        << endl;


    cout
        << "Total input rows:                  "
        << total_rows
        << endl;

    cout
        << "Peptide Category == protein:       "
        << protein_rows
        << endl;

    cout
        << "Pass length 50-2000 aa:            "
        << protein_length_pass
        << endl;

    cout
        << "Pass pLDDT >= 80:                  "
        << final_pass
        << endl;


    cout
        << endl;


    if (final_pass > 0) {

        cout
            << fixed
            << setprecision(2);

        cout
            << "Mean protein length:               "
            << mean_length
            << " aa"
            << endl;

        cout
            << "Protein length range:              "
            << min_selected_length
            << " - "
            << max_selected_length
            << " aa"
            << endl;

        cout
            << "Mean ColabFold pLDDT:              "
            << mean_plddt
            << endl;

        cout
            << "ColabFold pLDDT range:             "
            << min_selected_plddt
            << " - "
            << max_selected_plddt
            << endl;

        if (
            has_ptm &&
            valid_ptm_after_filter > 0
        ) {

            cout
                << "Mean ColabFold pTM:                "
                << mean_ptm
                << endl;
        }
    }


    cout
        << endl
        << "Output CSV:"
        << endl
        << output_csv
        << endl;


    cout
        << endl
        << "Summary:"
        << endl
        << summary_file
        << endl;


    return 0;
}