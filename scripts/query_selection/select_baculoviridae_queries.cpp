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
#include <map>

using namespace std;
namespace fs = std::filesystem;


// ============================================================
// trim
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
// lowercase
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
// ============================================================
string collapseSpaces(const string& s) {

    string result;
    bool previous_space = false;

    for (unsigned char c : s) {

        if (isspace(c)) {

            if (!previous_space) {
                result += ' ';
            }

            previous_space = true;

        } else {

            result += static_cast<char>(c);
            previous_space = false;
        }
    }

    return trim(result);
}


// ============================================================
// 标准化 Viro3D Name
// ============================================================
string normalizeName(string name) {

    name = trim(name);
    name = toLower(name);
    name = collapseSpaces(name);

    const string prefix = "product:";

    if (
        name.size() >= prefix.size() &&
        name.compare(0, prefix.size(), prefix) == 0
    ) {

        name = trim(
            name.substr(prefix.size())
        );
    }

    // 如果后面有 Note 等注释，只取主要 Product 名称
    size_t semicolon_pos =
        name.find(';');

    if (semicolon_pos != string::npos) {

        name =
            name.substr(
                0,
                semicolon_pos
            );
    }

    return trim(name);
}


// ============================================================
// TSV parser
// ============================================================
vector<string> parseTSVLine(
    const string& line
) {

    vector<string> fields;

    size_t start = 0;

    while (true) {

        size_t pos =
            line.find('\t', start);

        if (pos == string::npos) {

            fields.push_back(
                line.substr(start)
            );

            break;
        }

        fields.push_back(
            line.substr(
                start,
                pos - start
            )
        );

        start =
            pos + 1;
    }

    return fields;
}


// ============================================================
// 安全读取double
// ============================================================
bool parseDouble(
    const string& text,
    double& value
) {

    string x =
        trim(text);

    if (x.empty() ||
        x == "NA") {

        return false;
    }

    try {

        size_t pos = 0;

        value =
            stod(x, &pos);

        return
            pos == x.size();

    } catch (...) {

        return false;
    }
}


// ============================================================
// 一条候选query
// ============================================================
struct Candidate {

    vector<string> fields;

    string viro3d_id;
    string species;
    string raw_name;
    string normalized_name;
    string pdb_path;

    double plddt = 0.0;
    double ptm = 0.0;
};


// ============================================================
// 主程序
// ============================================================
int main() {

    // ========================================================
    // 输入
    // ========================================================

    const string input_file =
        "Baculoviridae_structure_mapping.tsv";


    // ========================================================
    // 输出目录
    // ========================================================

    const string output_dir =
        "processed/Baculoviridae/final";


    fs::create_directories(
        output_dir
    );


    const string selected_file =
        output_dir +
        "/selected_queries.tsv";


    const string id_file =
        output_dir +
        "/selected_query_ids.txt";


    const string pdb_file =
        output_dir +
        "/selected_query_pdb_paths.txt";


    const string summary_file =
        output_dir +
        "/query_selection_summary.txt";


    // ========================================================
    // 选定的12类
    // ========================================================

    const vector<string>
        target_functions = {

            "chitinase",
            "dna polymerase",
            "vp1054",
            "me53",
            "sod",
            "gp37",
            "ubiquitin",
            "lef-1",
            "lef-4",
            "p49",
            "p33",
            "egt"
        };


    const size_t
        QUERIES_PER_FUNCTION = 4;


    // ========================================================
    // 读取输入
    // ========================================================

    ifstream in(
        input_file
    );


    if (!in.is_open()) {

        cerr
            << "ERROR: cannot open:\n"
            << input_file
            << endl;

        return 1;
    }


    string line;


    if (!getline(in, line)) {

        cerr
            << "ERROR: empty input."
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
        parseTSVLine(
            line
        );


    unordered_map<string, size_t>
        col;


    for (
        size_t i = 0;
        i < headers.size();
        ++i
    ) {

        col[
            trim(headers[i])
        ] = i;
    }


    // ========================================================
    // 检查需要的列
    // ========================================================

    vector<string>
        required_columns = {

            "Viro3D_ID",
            "ICTV_Species",
            "Viro3D_Name",
            "ColabFold_pLDDT",
            "ColabFold_pTM",
            "Structure_Status",
            "Dataset_PDB_Path"
        };


    for (
        const string& name :
        required_columns
    ) {

        if (
            col.find(name)
            ==
            col.end()
        ) {

            cerr
                << "ERROR: missing column: "
                << name
                << endl;

            return 1;
        }
    }


    size_t idx_id =
        col["Viro3D_ID"];

    size_t idx_species =
        col["ICTV_Species"];

    size_t idx_name =
        col["Viro3D_Name"];

    size_t idx_plddt =
        col["ColabFold_pLDDT"];

    size_t idx_ptm =
        col["ColabFold_pTM"];

    size_t idx_status =
        col["Structure_Status"];

    size_t idx_pdb =
        col["Dataset_PDB_Path"];


    // ========================================================
    // function -> candidates
    // ========================================================

    unordered_map<string, vector<Candidate>>
        candidates;


    set<string>
        allowed_functions(
            target_functions.begin(),
            target_functions.end()
        );


    // ========================================================
    // 遍历所有记录
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
            parseTSVLine(
                line
            );


        if (
            fields.size()
            <
            headers.size()
        ) {

            continue;
        }


        string status =
            trim(
                fields[idx_status]
            );


        // 必须有结构
        if (
            status != "MATCHED"
        ) {

            continue;
        }


        string raw_name =
            trim(
                fields[idx_name]
            );


        string normalized_name =
            normalizeName(
                raw_name
            );


        // 只留下目标12类
        if (
            allowed_functions.find(
                normalized_name
            )
            ==
            allowed_functions.end()
        ) {

            continue;
        }


        double plddt = 0.0;
        double ptm = 0.0;


        if (
            !parseDouble(
                fields[idx_plddt],
                plddt
            )
        ) {

            continue;
        }


        if (
            !parseDouble(
                fields[idx_ptm],
                ptm
            )
        ) {

            continue;
        }


        Candidate c;

        c.fields =
            fields;

        c.viro3d_id =
            trim(
                fields[idx_id]
            );

        c.species =
            trim(
                fields[idx_species]
            );

        c.raw_name =
            raw_name;

        c.normalized_name =
            normalized_name;

        c.plddt =
            plddt;

        c.ptm =
            ptm;

        c.pdb_path =
            trim(
                fields[idx_pdb]
            );


        candidates[
            normalized_name
        ].push_back(
            c
        );
    }


    in.close();


    // ========================================================
    // 每一类选择4条
    //
    // 排序：
    // 1. pLDDT 高
    // 2. pTM 高
    //
    // 但是必须来自不同 species
    // ========================================================

    vector<Candidate>
        selected;


    map<string, vector<Candidate>>
        selected_by_function;


    for (
        const string& function :
        target_functions
    ) {

        auto it =
            candidates.find(
                function
            );


        if (
            it ==
            candidates.end()
        ) {

            cerr
                << "WARNING: no candidates for "
                << function
                << endl;

            continue;
        }


        vector<Candidate>& vec =
            it->second;


        sort(
            vec.begin(),
            vec.end(),

            [](
                const Candidate& a,
                const Candidate& b
            ) {

                if (
                    a.plddt !=
                    b.plddt
                ) {

                    return
                        a.plddt >
                        b.plddt;
                }


                if (
                    a.ptm !=
                    b.ptm
                ) {

                    return
                        a.ptm >
                        b.ptm;
                }


                return
                    a.viro3d_id <
                    b.viro3d_id;
            }
        );


        set<string>
            used_species;


        for (
            const Candidate& c :
            vec
        ) {

            if (
                used_species.find(
                    c.species
                )
                !=
                used_species.end()
            ) {

                continue;
            }


            selected.push_back(
                c
            );


            selected_by_function[
                function
            ].push_back(
                c
            );


            used_species.insert(
                c.species
            );


            if (
                used_species.size()
                >=
                QUERIES_PER_FUNCTION
            ) {

                break;
            }
        }
    }


    // ========================================================
    // 写 selected_queries.tsv
    // ========================================================

    ofstream selected_out(
        selected_file
    );


    // 原始metadata全部保留
    for (
        size_t i = 0;
        i < headers.size();
        ++i
    ) {

        if (i > 0) {
            selected_out << '\t';
        }

        selected_out
            << headers[i];
    }


    // 再增加标准化名称和选择顺序
    selected_out
        << '\t'
        << "Standardized_Function"
        << '\t'
        << "Query_Group_Rank"
        << '\n';


    map<string, size_t>
        group_rank;


    for (
        const Candidate& c :
        selected
    ) {

        for (
            size_t i = 0;
            i < c.fields.size();
            ++i
        ) {

            if (i > 0) {
                selected_out << '\t';
            }

            selected_out
                << c.fields[i];
        }


        group_rank[
            c.normalized_name
        ]++;


        selected_out
            << '\t'
            << c.normalized_name
            << '\t'
            << group_rank[
                c.normalized_name
            ]
            << '\n';
    }


    selected_out.close();


    // ========================================================
    // 输出 query IDs
    // ========================================================

    ofstream id_out(
        id_file
    );


    for (
        const Candidate& c :
        selected
    ) {

        id_out
            << c.viro3d_id
            << '\n';
    }


    id_out.close();


    // ========================================================
    // 输出 PDB paths
    // ========================================================

    ofstream pdb_out(
        pdb_file
    );


    for (
        const Candidate& c :
        selected
    ) {

        pdb_out
            << c.pdb_path
            << '\n';
    }


    pdb_out.close();


    // ========================================================
    // summary
    // ========================================================

    ofstream summary(
        summary_file
    );


    summary
        << "Baculoviridae query selection summary\n";

    summary
        << "========================================\n\n";


    summary
        << "Reference proteins: 6843\n";

    summary
        << "Target functional groups: "
        << target_functions.size()
        << "\n";

    summary
        << "Queries per functional group: "
        << QUERIES_PER_FUNCTION
        << "\n";

    summary
        << "Total selected queries: "
        << selected.size()
        << "\n\n";


    summary
        << "Selection rules:\n";

    summary
        << "  1. Structure_Status == MATCHED\n";

    summary
        << "  2. Select from 12 predefined functional groups\n";

    summary
        << "  3. Four proteins per functional group\n";

    summary
        << "  4. Proteins within the same group must come from different ICTV species\n";

    summary
        << "  5. Higher ColabFold pLDDT is preferred\n";

    summary
        << "  6. Higher ColabFold pTM is used as the secondary ranking criterion\n\n";


    summary
        << "Selected groups:\n";


    for (
        const string& function :
        target_functions
    ) {

        summary
            << "\n"
            << function
            << "\n";


        const auto& group =
            selected_by_function[
                function
            ];


        for (
            const Candidate& c :
            group
        ) {

            summary
                << "  "
                << c.viro3d_id
                << "\t"
                << c.species
                << "\tpLDDT="
                << fixed
                << setprecision(2)
                << c.plddt
                << "\tpTM="
                << c.ptm
                << "\n";
        }
    }


    summary.close();


    // ========================================================
    // 打印结果
    // ========================================================

    cout
        << "========================================"
        << endl;

    cout
        << "Query selection"
        << endl;

    cout
        << "========================================"
        << endl;


    cout
        << "Functional groups: "
        << target_functions.size()
        << endl;

    cout
        << "Queries per group: "
        << QUERIES_PER_FUNCTION
        << endl;

    cout
        << "Total queries: "
        << selected.size()
        << endl;


    cout
        << endl;


    for (
        const string& function :
        target_functions
    ) {

        cout
            << function
            << ": "
            << selected_by_function[
                function
            ].size()
            << " queries"
            << endl;
    }


    cout
        << endl
        << "Selected queries:"
        << endl
        << selected_file
        << endl;


    cout
        << endl
        << "Query IDs:"
        << endl
        << id_file
        << endl;


    cout
        << endl
        << "Query PDB paths:"
        << endl
        << pdb_file
        << endl;


    cout
        << endl
        << "Summary:"
        << endl
        << summary_file
        << endl;


    return 0;
}