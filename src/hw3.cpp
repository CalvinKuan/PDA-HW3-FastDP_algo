#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <sstream>
#include <cctype>

using namespace std;
using Clock = std::chrono::high_resolution_clock;
using ns = std::chrono::nanoseconds;

// ---- 全域 timing 累計（nanoseconds）----
static long long g_t_opt_region_ns = 0;  // compute_optimal_region_for_inst
static long long g_t_legal_same_ns = 0;  // check_legal_swap_same_row
static long long g_t_legal_two_ns = 0;   // check_legal_swap_two_rows
static long long g_t_delta_hpwl_ns = 0;  // delta_hpwl_cached
static long long g_t_hpwl_recalc_ns = 0; // calculate_HPWL（兩個 for 迴圈）

// ---- 呼叫次數（方便算平均）----
static long long g_cnt_opt_region = 0;
static long long g_cnt_legal_same = 0;
static long long g_cnt_legal_two = 0;
static long long g_cnt_delta_hpwl = 0;
static long long g_cnt_hpwl_recalc = 0;

// （可選）重設 timing 的小工具
static void reset_global_swap_timers()
{
    g_t_opt_region_ns = 0;
    g_t_legal_same_ns = 0;
    g_t_legal_two_ns = 0;
    g_t_delta_hpwl_ns = 0;
    g_t_hpwl_recalc_ns = 0;

    g_cnt_opt_region = 0;
    g_cnt_legal_same = 0;
    g_cnt_legal_two = 0;
    g_cnt_delta_hpwl = 0;
    g_cnt_hpwl_recalc = 0;
}

// （可選）印出統計結果
static void print_global_swap_timers()
{
    auto to_ms = [](long long ns_val)
    {
        return ns_val / 1e6; // 轉毫秒
    };

    std::fprintf(stderr,
                 "[GS timing] opt_region  : %.3f ms (%lld calls)\n"
                 "[GS timing] legal_same  : %.3f ms (%lld calls)\n"
                 "[GS timing] legal_two   : %.3f ms (%lld calls)\n"
                 "[GS timing] delta_hpwl  : %.3f ms (%lld calls)\n"
                 "[GS timing] hpwl_recalc : %.3f ms (%lld calls)\n",
                 to_ms(g_t_opt_region_ns), g_cnt_opt_region,
                 to_ms(g_t_legal_same_ns), g_cnt_legal_same,
                 to_ms(g_t_legal_two_ns), g_cnt_legal_two,
                 to_ms(g_t_delta_hpwl_ns), g_cnt_delta_hpwl,
                 to_ms(g_t_hpwl_recalc_ns), g_cnt_hpwl_recalc);
}

int dbu;
int site_width_dbu;
int TIME_LIMIT_MS;

struct macro_info
{
    string class_type;
    double size_x_micron;
    double size_y_micron;
    string site;
};
struct row_info
{
    string row_name;
    string site;
    int origin_x_dbu;
    int origin_y_dbu;
    string orient;
    int x_site_size;
    int y_site_size;
    int x_spacing_dbu;
};
struct inst_info
{
    string id;
    string type;
    string mov_or_fix;
    int x_dbu;
    int y_dbu;
    string orient;
};
struct pin_map
{
    string net_name;
    int x_dbu;
    int y_dbu;
    string orient;
};
struct net_map
{
    vector<string> connected_insts;
};
struct RowStripe
{
    string row_name;
    int x_origin_dbu;
    int y_origin_dbu;
    int x_sites;
    string orient;
    int x_spacing_dbu;
    vector<pair<int, int>> row_blocks; // [xL,xR) intervals blocked by macros
    vector<inst_info> placed_insts;
};
struct rect
{
    int x_left;
    int x_right;
    int y_bottom;
    int y_top;
};

void print_row_with_cells(vector<RowStripe> &row_stripes)
{
    int count = 0;
    for (auto &row : row_stripes)
    {

        cout << "Row:" << row.row_name << ": Origin(" << row.x_origin_dbu << ", " << row.y_origin_dbu << "), Sites: " << row.x_sites;
        cout << "  Placed Insts: " << endl;
        for (auto &inst : row.placed_insts)
        {
            count++;
            cout << inst.id << " " << "(" << inst.x_dbu << ", " << inst.y_dbu << ") " << endl;
        }
        cout << "  Blocked Regions: ";
        for (auto &block : row.row_blocks)
        {
            cout << "(" << block.first << ", " << block.second << ") ";
        }
        cout << endl;

        cout << endl;
    }
    // cout << "count:" <<count<<endl;
}

//---------------- helpers ----------------//
vector<RowStripe> make_cell_map(vector<row_info> &row_table,
                                unordered_map<string, inst_info> &inst_table,
                                unordered_map<string, inst_info> &block_map)
{
    vector<RowStripe> row_stripes;
    row_stripes.resize(row_table.size());

    // 先把 row_table 的資訊塞進 row_stripes，並建 y 座標 -> row index 的 map
    unordered_map<int, int> y2row;
    y2row.reserve(row_table.size() * 2);

    for (int i = 0; i < (int)row_table.size(); ++i)
    {
        auto &row = row_table[i];
        RowStripe rs;
        rs.row_name = row.row_name;
        rs.x_origin_dbu = row.origin_x_dbu;
        rs.y_origin_dbu = row.origin_y_dbu;
        rs.x_sites = row.x_site_size;
        rs.x_spacing_dbu = row.x_spacing_dbu;
        rs.orient = row.orient;
        rs.row_blocks.clear();
        rs.placed_insts.clear();

        row_stripes[i] = std::move(rs);
        y2row[row.origin_y_dbu] = i;
    }

    for (auto &kv : inst_table)
    {
        const inst_info &inst = kv.second;
        if (inst.mov_or_fix != "PLACED")
            continue;

        auto it = y2row.find(inst.y_dbu);
        if (it == y2row.end())
            continue; // 找不到對應 row 的就跳過（理論上不會發生）

        int ridx = it->second;
        row_stripes[ridx].placed_insts.push_back(inst);
    }

    // 每一 row 按 x 排序一次
    for (auto &row : row_stripes)
    {
        std::sort(row.placed_insts.begin(), row.placed_insts.end(),
                  [](const inst_info &a, const inst_info &b)
                  { return a.x_dbu < b.x_dbu; });
    }

    return row_stripes;
}

void makey2row_idx(vector<RowStripe> &rowstripes, unordered_map<int, int> &y2row_idx)
{
    for (int i = 0; i < rowstripes.size(); i++)
    {
        y2row_idx.insert({rowstripes[i].y_origin_dbu, i});
    }
}

void calculate_block_bound(unordered_map<string, inst_info> &block_map,
                           unordered_map<string, macro_info> &library,
                           vector<RowStripe> &row_stripes)
{
    for (auto &kv : block_map)
    {
        const inst_info &blk = kv.second;
        auto it = library.find(blk.type);
        if (it == library.end())
            continue;
        const auto &mac = it->second;
        int xL = blk.x_dbu;
        int xR = blk.x_dbu + (int)llround(mac.size_x_micron * dbu);
        int yB = blk.y_dbu;
        int yT = blk.y_dbu + (int)llround(mac.size_y_micron * dbu);
        for (auto &row : row_stripes)
        {
            if (row.y_origin_dbu >= yB && row.y_origin_dbu < yT)
            {
                row.row_blocks.push_back({xL, xR});
            }
        }
    }
}

static long long delta_hpwl_cached(
    const string &A, const string &B,
    unordered_map<string, inst_info> &inst_table,
    unordered_map<string, net_map> &net_table,
    unordered_map<string, vector<string>> &inst2nets,
    unordered_map<string, pin_map> &pin_table,
    unordered_map<string, long long> &net2hpwl,
    int ax, int ay, int bx, int by,
    std::vector<std::pair<std::string, long long>> *out_new = nullptr)
{
    long long delta = 0;

    // 用 static 小 vector 重用 buffer，避免每次 allocate/free
    static std::vector<std::string> nets;
    nets.clear();

    auto &na = inst2nets[A];
    auto &nb = inst2nets[B];
    nets.reserve(na.size() + nb.size());
    nets.insert(nets.end(), na.begin(), na.end());
    nets.insert(nets.end(), nb.begin(), nb.end());

    // 去重：資料量通常很小，用 sort + unique 反而比 unordered_set 快
    std::sort(nets.begin(), nets.end());
    nets.erase(std::unique(nets.begin(), nets.end()), nets.end());

    if (out_new)
        out_new->clear();

    for (const auto &netname : nets)
    {
        long long before = net2hpwl[netname];

        const auto &net = net_table[netname];
        int x_max = INT_MIN, x_min = INT_MAX;
        int y_max = INT_MIN, y_min = INT_MAX;

        for (const auto &name : net.connected_insts)
        {
            int xx, yy;
            if (name == A)
            {
                xx = ax;
                yy = ay;
            }
            else if (name == B)
            {
                xx = bx;
                yy = by;
            }
            else
            {
                auto itI = inst_table.find(name);
                if (itI != inst_table.end())
                {
                    xx = itI->second.x_dbu;
                    yy = itI->second.y_dbu;
                }
                else
                {
                    const auto &p = pin_table[name];
                    xx = p.x_dbu;
                    yy = p.y_dbu;
                }
            }
            x_max = max(x_max, xx);
            x_min = min(x_min, xx);
            y_max = max(y_max, yy);
            y_min = min(y_min, yy);
        }

        long long after = (x_max - x_min) + (y_max - y_min);
        delta += (after - before);

        if (out_new)
            out_new->emplace_back(netname, after);
    }

    return delta;
}


// HPWL for one net list (instances or IO pins)
long long calculate_HPWL(const vector<string> &net,
                         unordered_map<string, inst_info> &inst_table,
                         unordered_map<string, pin_map> &pin_table)
{
    long long hpwl = 0;
    int min_x = INT_MAX;
    int max_x = INT_MIN;
    int min_y = INT_MAX;
    int max_y = INT_MIN;
    bool any = false;
    for (const auto &name : net)
    {
        auto itI = inst_table.find(name);
        if (itI != inst_table.end())
        {

            const inst_info &inst = itI->second;

            min_x = min(min_x, inst.x_dbu);
            max_x = max(max_x, inst.x_dbu);
            min_y = min(min_y, inst.y_dbu);
            max_y = max(max_y, inst.y_dbu);
            any = true;
        }
        else
        {
            auto itP = pin_table.find(name);
            if (itP != pin_table.end())
            {
                const pin_map &pin = itP->second;
                min_x = min(min_x, pin.x_dbu);
                max_x = max(max_x, pin.x_dbu);
                min_y = min(min_y, pin.y_dbu);
                max_y = max(max_y, pin.y_dbu);
                any = true;
            }
        }
    }
    if (!any)
        return 0;
    hpwl = (long long)(max_x - min_x) + (long long)(max_y - min_y);
    return hpwl;
}

static unordered_map<string, long long> build_net2hpwl_cache(
    unordered_map<string, net_map> &net_table,
    unordered_map<string, inst_info> &inst_table,
    unordered_map<string, pin_map> &pin_table)
{
    unordered_map<string, long long> net2hpwl;
    int count = 0;
    for (auto &kv : net_table)
    {
        // cout<< kv.first <<" "<<net2hpwl[kv.first]<<endl;
        net2hpwl[kv.first] = calculate_HPWL(kv.second.connected_insts, inst_table, pin_table);
        // cout<< kv.first <<" "<<net2hpwl[kv.first]<<endl;
        count++;
    }
    // cout<<count<<endl;
    return net2hpwl;
}

long long calculate_total_HPWL(unordered_map<string, net_map> &net_table,
                               unordered_map<string, inst_info> &inst_table,
                               unordered_map<string, pin_map> &pin_table)
{
    long long total_hpwl = 0;
    for (auto &kv : net_table)
    {
        total_hpwl += calculate_HPWL(kv.second.connected_insts, inst_table, pin_table);
    }
    return total_hpwl;
}

unordered_map<string, vector<string>> build_inst2nets(const unordered_map<string, net_map> &net_table)
{
    unordered_map<string, vector<string>> m;
    for (auto &kv : net_table)
    {
        const string &netname = kv.first;
        for (auto &name : kv.second.connected_insts)
            m[name].push_back(netname);
    }
    return m;
}

unordered_map<string, int> build_macro_width_dbu(const unordered_map<string, macro_info> &lib)
{
    unordered_map<string, int> w;
    for (auto &kv : lib)
    {
        const auto &mac = kv.second;
        int wx = (int)llround(mac.size_x_micron * dbu);
        w[kv.first] = wx;
    }
    return w;
}

// compute BBox excluding "inst" itself
static bool bbox_excluding(const inst_info &me,
                           unordered_map<string, net_map> &net_table,
                           unordered_map<string, pin_map> &pin_table,
                           const string &net_name,
                           unordered_map<string, inst_info> &inst_table,
                           vector<int> &L, vector<int> &R, vector<int> &T, vector<int> &B)
{
    auto it = net_table.find(net_name);
    if (it == net_table.end())
        return false;
    const auto &net = it->second;
    int x_left = INT_MAX, x_right = INT_MIN, y_bottom = INT_MAX, y_top = INT_MIN;
    bool flag = false;
    for (const auto &name : net.connected_insts)
    {
        if (name == me.id)
            continue;
        auto itI = inst_table.find(name);
        if (itI != inst_table.end())
        {
            const auto &o = itI->second;
            x_left = min(x_left, o.x_dbu);
            x_right = max(x_right, o.x_dbu);
            y_bottom = min(y_bottom, o.y_dbu);
            y_top = max(y_top, o.y_dbu);
            flag = true;
            continue;
        }
        auto itP = pin_table.find(name);
        if (itP != pin_table.end())
        {
            const auto &p = itP->second;
            x_left = min(x_left, p.x_dbu);
            x_right = max(x_right, p.x_dbu);
            y_bottom = min(y_bottom, p.y_dbu);
            y_top = max(y_top, p.y_dbu);
            flag = true;
        }
    }
    if (!flag)
        return false;
    L.push_back(x_left);
    R.push_back(x_right);
    B.push_back(y_bottom);
    T.push_back(y_top);
    return true;
}

static int median_value(vector<int> &v)
{
    if (v.empty())
        return 0;
    size_t mid = v.size() / 2;
    nth_element(v.begin(), v.begin() + mid, v.end());
    return v[mid];
}

static rect cal_optimal_region(vector<int> &L, vector<int> &R, vector<int> &T, vector<int> &B)
{
    rect optimal_region;
    int iL = median_value(L);
    int iR = median_value(R);
    int iB = median_value(B);
    int iT = median_value(T);
    optimal_region.x_left = min(iL, iR);
    optimal_region.x_right = max(iL, iR);
    optimal_region.y_bottom = min(iB, iT);
    optimal_region.y_top = max(iB, iT);

    // pad region to avoid over-narrow boxes
    optimal_region.x_left;
    optimal_region.x_right;
    optimal_region.y_bottom;
    optimal_region.y_top;
    return optimal_region;
}

static rect compute_optimal_region_for_inst(const string &inst_name,
                                            unordered_map<string, net_map> &net_table,
                                            unordered_map<string, pin_map> &pin_table,
                                            unordered_map<string, inst_info> &inst_table,
                                            unordered_map<string, vector<string>> &inst2nets)
{
    vector<int> L, R, T, B;
    auto itc = inst_table.find(inst_name);
    const inst_info &me = itc->second;
    auto itn = inst2nets.find(inst_name);
    if (itn != inst2nets.end())
    {
        for (const auto &netname : itn->second)
        {
            bbox_excluding(me, net_table, pin_table, netname, inst_table, L, R, T, B);
        }
    }
    return cal_optimal_region(L, R, T, B);
}

//---------------- legality helpers ----------------//
static bool no_overlap_with_blocks(int L, int R, const vector<pair<int, int>> &blocks)
{
    for (const auto &b : blocks)
    {
        if (!(R <= b.first || b.second <= L))
            return false;
    }
    return true;
}
static bool inside_row_usable_range(int L, int R, int row_xL, int row_xR)
{
    return (L >= row_xL) && (R <= row_xR);
}
static bool is_site_aligned(int x, int row_xL, int site_w)
{
    if (site_w <= 0)
        return true;
    return ((x - row_xL) % site_w) == 0;
}

// full-scan overlap check on a row
static bool no_overlap_with_row(const RowStripe &row,
                                const unordered_map<string, int> &width_dbu,
                                int newL, int newR,
                                int idx_ignore1, int idx_ignore2)
{
    for (int i = 0; i < (int)row.placed_insts.size(); ++i)
    {
        if (i == idx_ignore1 || i == idx_ignore2)
            continue;
        const auto &z = row.placed_insts[i];
        auto itW = width_dbu.find(z.type);
        if (itW == width_dbu.end())
            return false;
        int zL = z.x_dbu, zR = z.x_dbu + itW->second;
        if (!(newR <= zL || zR <= newL))
            return false;
    }
    return true;
}

static bool check_legal_swap_same_row(
    const RowStripe &rowstripes,
    int ia, int ib,
    const unordered_map<string, int> &width_dbu,
    int row_xL, int row_xR,
    int site_w)
{
    const inst_info &A = rowstripes.placed_insts[ia];
    const inst_info &B = rowstripes.placed_insts[ib];

    auto itA = width_dbu.find(A.type);
    auto itB = width_dbu.find(B.type);
    if (itA == width_dbu.end() || itB == width_dbu.end())
        return false;
    int wA = itA->second;
    int wB = itB->second;

    int aL = B.x_dbu, aR = B.x_dbu + wA;
    int bL = A.x_dbu, bR = A.x_dbu + wB;

    if (!no_overlap_with_blocks(aL, aR, rowstripes.row_blocks))
        return false;
    if (!no_overlap_with_blocks(bL, bR, rowstripes.row_blocks))
        return false;
    if (!inside_row_usable_range(aL, aR, row_xL, row_xR))
        return false;
    if (!inside_row_usable_range(bL, bR, row_xL, row_xR))
        return false;
    if (!is_site_aligned(aL, row_xL, site_w))
        return false;
    if (!is_site_aligned(bL, row_xL, site_w))
        return false;

    if (!no_overlap_with_row(rowstripes, width_dbu, aL, aR, ia, ib))
        return false;
    if (!no_overlap_with_row(rowstripes, width_dbu, bL, bR, ia, ib))
        return false;

    if (abs(ia - ib) == 1)
    {
        if (ia < ib)
        {
            if (bR > aL)
            {
                return false;
            }
        }
        else if (ib < ia)
        {
            if (aR > bL)
            {
                return false;
            }
        }
    }
    return true;
}

static bool check_legal_swap_two_rows(
    const RowStripe &rowA, int ia, int rowA_xL, int rowA_xR,
    const RowStripe &rowB, int ib, int rowB_xL, int rowB_xR,
    const unordered_map<string, int> &width_dbu,
    int site_w)
{
    const auto &A = rowA.placed_insts[ia];
    const auto &B = rowB.placed_insts[ib];

    auto itA = width_dbu.find(A.type);
    auto itB = width_dbu.find(B.type);
    if (itA == width_dbu.end() || itB == width_dbu.end())
        return false;
    int wA = itA->second, wB = itB->second;

    int aL = B.x_dbu, aR = B.x_dbu + wA; // A -> rowB
    int bL = A.x_dbu, bR = A.x_dbu + wB; // B -> rowA

    if (!no_overlap_with_blocks(aL, aR, rowB.row_blocks))
        return false;
    if (!no_overlap_with_blocks(bL, bR, rowA.row_blocks))
        return false;
    if (!inside_row_usable_range(aL, aR, rowB_xL, rowB_xR))
        return false;
    if (!inside_row_usable_range(bL, bR, rowA_xL, rowA_xR))
        return false;
    if (!is_site_aligned(aL, rowB_xL, site_w))
        return false;
    if (!is_site_aligned(bL, rowA_xL, site_w))
        return false;

    if (!no_overlap_with_row(rowB, width_dbu, aL, aR, ib, -1))
        return false;
    if (!no_overlap_with_row(rowA, width_dbu, bL, bR, ia, -1))
        return false;

    return true;
}

// 檢查：在 row 的 idx 這個位置放一顆 new cell (newL,newR)
// 會不會超出 row 邊界 / 撞到 block / 撞到左右鄰居
static bool ok_row_single_swap_with_blocks(
    const RowStripe &row,
    int idx,                                                  // 要放新 cell 的 index（原本就有 cell 在這裡）
    int newL, int newR, unordered_map<string, int> &width_dbu // 新 cell 的 [xL, xR)
)
{
    // ---- 1) row 左右邊界 ----
    int rowL = row.x_origin_dbu;
    int rowR = row.x_origin_dbu + row.x_sites * row.x_spacing_dbu;
    if (newL < rowL || newR > rowR)
        return false;

    // ---- 2) block 區段檢查 ----
    // row_blocks 裡是 [blkL, blkR) 的 blocked 區域
    for (const auto &blk : row.row_blocks)
    {
        int blkL = blk.first;
        int blkR = blk.second;
        // 只要新區間與 block 有交集就不行
        if (!(newR <= blkL || blkR <= newL))
            return false;
    }

    // ---- 3) 左右鄰居檢查（placed_insts 已經是 legal / 排序）----
    const int N = (int)row.placed_insts.size();

    // 左鄰居
    if (idx - 1 >= 0)
    {
        const inst_info &L = row.placed_insts[idx - 1];
        int Lw = width_dbu[L.type]; // 如果你沒有這個欄位，改用 width_dbu[L.type]
        int Ll = L.x_dbu;
        int Lr = Ll + Lw;

        // new 區間左端必須在左鄰居右邊
        if (!(newL >= Lr))
            return false;
    }

    // 右鄰居
    if (idx + 1 < N)
    {
        const inst_info &R = row.placed_insts[idx + 1];
        int Rw = width_dbu[R.type]; // 同上，如果沒有就 width_dbu[R.type]
        int Rl = R.x_dbu;
        int Rr = Rl + Rw;

        // new 區間右端必須在右鄰居左邊
        if (!(newR <= Rl))
            return false;
    }

    return true;
}

static bool check_legal_swap_two_rows_fast(
    const RowStripe &rowA, int idxA,
    const RowStripe &rowB, int idxB,
    unordered_map<string, int> &width_dbu)
{
    const inst_info &A = rowA.placed_insts[idxA];
    const inst_info &B = rowB.placed_insts[idxB];

    // 取得寬度（如果 inst_info 本身有 width，就可以不用查 map）
    int wA, wB;
    auto itWA = width_dbu.find(A.type);
    auto itWB = width_dbu.find(B.type);
    if (itWA == width_dbu.end() || itWB == width_dbu.end())
        return false;
    wA = itWA->second;
    wB = itWB->second;

    // A 現在在 rowA 的 x = A.x_dbu，B 在 rowB 的 x = B.x_dbu
    // swap 後：
    //  - rowA 的 idxA 位置放 B，區間 [A.x, A.x + wB)
    //  - rowB 的 idxB 位置放 A，區間 [B.x, B.x + wA)

    int newBL = A.x_dbu;
    int newBR = newBL + wB;

    int newAL = B.x_dbu;
    int newAR = newAL + wA;

    // rowA 上放 B 會不會超界 / 撞 block / 撞左右鄰居
    if (!ok_row_single_swap_with_blocks(rowA, idxA, newBL, newBR, width_dbu))
        return false;

    // rowB 上放 A 會不會超界 / 撞 block / 撞左右鄰居
    if (!ok_row_single_swap_with_blocks(rowB, idxB, newAL, newAR, width_dbu))
        return false;

    return true;
}

// delta HPWL if A placed at (ax,ay) and B placed at (bx,by)
static long long delta_hpwl(const string &inst_a_name, const string &inst_b_name,
                            unordered_map<string, inst_info> &inst_table,
                            unordered_map<string, net_map> &net_table,
                            unordered_map<string, vector<string>> &inst2nets,
                            unordered_map<string, pin_map> &pin_table,
                            int ax, int ay, int bx, int by)
{
    long long delta = 0;
    unordered_set<string> affected_nets;
    auto &nets_with_a = inst2nets[inst_a_name];
    auto &nets_with_b = inst2nets[inst_b_name];
    affected_nets.insert(nets_with_a.begin(), nets_with_a.end());
    affected_nets.insert(nets_with_b.begin(), nets_with_b.end());

    for (const auto &netname : affected_nets)
    {
        auto &net = net_table[netname];
        int before_swap = (int)calculate_HPWL(net.connected_insts, inst_table, pin_table);

        int x_max = INT_MIN, x_min = INT_MAX, y_max = INT_MIN, y_min = INT_MAX;
        for (const auto &name : net.connected_insts)
        {
            if (name == inst_a_name)
            {
                x_max = max(x_max, ax);
                x_min = min(x_min, ax);
                y_max = max(y_max, ay);
                y_min = min(y_min, ay);
            }
            else if (name == inst_b_name)
            {
                x_max = max(x_max, bx);
                x_min = min(x_min, bx);
                y_max = max(y_max, by);
                y_min = min(y_min, by);
            }
            else
            {
                auto itI = inst_table.find(name);
                if (itI != inst_table.end())
                {
                    const auto &inst = itI->second;
                    x_max = max(x_max, inst.x_dbu);
                    x_min = min(x_min, inst.x_dbu);
                    y_max = max(y_max, inst.y_dbu);
                    y_min = min(y_min, inst.y_dbu);
                }
                else
                {
                    auto itP = pin_table.find(name);
                    if (itP != pin_table.end())
                    {
                        const auto &pin = itP->second;
                        x_max = max(x_max, pin.x_dbu);
                        x_min = min(x_min, pin.x_dbu);
                        y_max = max(y_max, pin.y_dbu);
                        y_min = min(y_min, pin.y_dbu);
                    }
                }
            }
        }
        int after_swap = (x_max - x_min) + (y_max - y_min);
        delta += (after_swap - before_swap);
    }
    return delta;
}

//---------------- Global Swap core ----------------//
static bool global_swap_one_pass(vector<RowStripe> &rowstripes,
                                 unordered_map<string, inst_info> &inst_table,
                                 unordered_map<string, net_map> &net_table,
                                 unordered_map<string, vector<string>> &inst2nets,
                                 unordered_map<string, pin_map> &pin_table,
                                 unordered_map<string, int> &width_dbu,
                                 int max_row_delta,
                                 long long &total_hpwl,
                                 unordered_map<string, long long> &net2hpwl,
                                 int K,
                                 int x_w,
                                 unordered_map<int, int> &y2row_idx)
{
    bool any_improved = false;
    const int ROWS = (int)rowstripes.size();

    // localized search params
    const int MAX_ROW_DELTA = max_row_delta; // 限制 r 上下的 row 數
    const int K_NEAR = K;                    // neighbor 數
    int X_WIN;                               // A.x ± window

    // 依設計大小動態決定 stride（跳著選 A）
    int total_insts = (int)inst_table.size();
    int stride = 1;
    if (total_insts > 800000)
        stride = 4;
    else if (total_insts > 400000)
        stride = 3;
    else if (total_insts > 200000)
        stride = 2;

    for (int r = 0; r < ROWS; ++r)
    {
        auto &row = rowstripes[r];
        X_WIN = x_w * row.x_spacing_dbu;
        const int N = (int)row.placed_insts.size();
        if (N == 0)
            continue;

        // 這裡用 (r % stride) 讓不同 row 的起點錯開一點
        for (int ia = (r % stride); ia < N; ia += stride)
        {
            inst_info &A = row.placed_insts[ia];

            // --- 計時 compute_optimal_region_for_inst ---
            auto t_opt0 = Clock::now();
            rect region = compute_optimal_region_for_inst(
                A.id, net_table, pin_table, inst_table, inst2nets);
            auto t_opt1 = Clock::now();
            g_t_opt_region_ns += std::chrono::duration_cast<ns>(t_opt1 - t_opt0).count();
            ++g_cnt_opt_region;
            // -------------------------------------------

            long long best_delta = 0;
            int best_r = -1, best_ib = -1;

            // 由 optimal region 查出大概的 row 範圍
            auto it_low = y2row_idx.find(region.y_bottom);
            auto it_high = y2row_idx.find(region.y_top);
            if (it_low == y2row_idx.end() || it_high == y2row_idx.end())
                continue;

            int rr_low = it_low->second;
            int rr_high = it_high->second;

            // 再用 MAX_ROW_DELTA 進一步收窄：只看 r ± MAX_ROW_DELTA
            int rr_start = std::max(rr_low, std::max(0, r - MAX_ROW_DELTA));
            int rr_end = std::min(rr_high, std::min(ROWS - 1, r + MAX_ROW_DELTA));

            for (int rr = rr_start; rr <= rr_end; ++rr)
            {
                auto &row2 = rowstripes[rr];
                auto &vec = row2.placed_insts;
                if (vec.empty())
                    continue;

                // region window（保險起見再檢查一次 y）
                if (row2.y_origin_dbu < region.y_bottom ||
                    row2.y_origin_dbu > region.y_top)
                    continue;

                // binary search around A.x
                auto it = std::lower_bound(
                    vec.begin(), vec.end(), A.x_dbu,
                    [](const inst_info &z, int x)
                    { return z.x_dbu < x; });
                int center = (int)std::distance(vec.begin(), it);
                int left = std::max(0, center - K_NEAR);
                int right = std::min((int)vec.size() - 1, center + K_NEAR);

                for (int ib = left; ib <= right; ++ib)
                {
                    const inst_info &B = vec[ib];
                    if (B.id == A.id)
                        continue;
                    if (std::abs(B.x_dbu - A.x_dbu) > X_WIN)
                        continue;

                    // region window (x / y)
                    if (B.x_dbu < region.x_left || B.x_dbu > region.x_right)
                        continue;
                    if (B.y_dbu < region.y_bottom || B.y_dbu > region.y_top)
                        continue;

                    bool legal = false;

                    auto t_legal0 = Clock::now();
                    if (r == rr)
                    {
                        // 你原本的 same-row 檢查（如果它已有處理 block 就照用）
                        legal = check_legal_swap_same_row(
                            row, ia, ib, width_dbu,
                            row.x_origin_dbu,
                            row.x_origin_dbu + row.x_sites * row.x_spacing_dbu,
                            row.x_spacing_dbu);
                        auto t_legal1 = Clock::now();
                        g_t_legal_same_ns += std::chrono::duration_cast<ns>(t_legal1 - t_legal0).count();
                        ++g_cnt_legal_same;
                    }
                    else
                    {
                        // 新的 two-row fast 版本（有看 row_blocks）
                        legal = check_legal_swap_two_rows_fast(
                            row, ia, row2, ib, width_dbu);
                        auto t_legal1 = Clock::now();
                        g_t_legal_two_ns += std::chrono::duration_cast<ns>(t_legal1 - t_legal0).count();
                        ++g_cnt_legal_two;
                    }

                    // -------------------------

                    if (!legal)
                        continue;

                    int ax = B.x_dbu;
                    int ay = (r == rr) ? A.y_dbu : row2.y_origin_dbu;
                    int bx = A.x_dbu;
                    int by = (r == rr) ? B.y_dbu : row.y_origin_dbu;

                    // --- 計時 delta_hpwl_cached ---
                    auto t_delta0 = Clock::now();
                    long long d = delta_hpwl_cached(
                        A.id, B.id,
                        inst_table, net_table, inst2nets, pin_table,
                        net2hpwl, ax, ay, bx, by);
                    auto t_delta1 = Clock::now();
                    g_t_delta_hpwl_ns += std::chrono::duration_cast<ns>(t_delta1 - t_delta0).count();
                    ++g_cnt_delta_hpwl;
                    // --------------------------------

                    if (d < best_delta)
                    {
                        best_delta = d;
                        best_r = rr;
                        best_ib = ib;
                    }
                }
            }

            if (best_delta < 0 && best_r >= 0)
            {
                any_improved = true;

                auto &rowB = rowstripes[best_r];
                auto &B = rowB.placed_insts[best_ib];

                int Ax = B.x_dbu, Ay = B.y_dbu;
                int Bx = A.x_dbu, By = A.y_dbu;

                // 更新 inst_table 的位置 + orient
                inst_table[A.id].x_dbu = Ax;
                inst_table[A.id].y_dbu = Ay;
                inst_table[B.id].x_dbu = Bx;
                inst_table[B.id].y_dbu = By;
                std::string tmp_orient = inst_table[A.id].orient;
                inst_table[A.id].orient = inst_table[B.id].orient;
                inst_table[B.id].orient = tmp_orient;

                inst_info a_new = inst_table[A.id];
                inst_info b_new = inst_table[B.id];

                // row vectors 也要同步
                row.placed_insts[ia] = b_new;
                if (best_r == r)
                    row.placed_insts[best_ib] = a_new;
                else
                    rowB.placed_insts[best_ib] = a_new;

                // 更新 total_hpwl
                total_hpwl += best_delta;

                // --- 計時 HPWL 重新計算（兩個 for 迴圈一起算） ---
                auto t_hpwl0 = Clock::now();

                for (const auto &netname : inst2nets[A.id])
                {
                    long long newv = calculate_HPWL(
                        net_table[netname].connected_insts,
                        inst_table, pin_table);
                    net2hpwl[netname] = newv;
                }
                for (const auto &netname : inst2nets[B.id])
                {
                    long long newv = calculate_HPWL(
                        net_table[netname].connected_insts,
                        inst_table, pin_table);
                    net2hpwl[netname] = newv;
                }

                auto t_hpwl1 = Clock::now();
                g_t_hpwl_recalc_ns += std::chrono::duration_cast<ns>(t_hpwl1 - t_hpwl0).count();
                ++g_cnt_hpwl_recalc;
                // ------------------------------------------------
            }
        }
    }

    return any_improved;
}

static bool vertical_swap_one_pass(
    vector<RowStripe> &rowstripes,
    unordered_map<string, inst_info> &inst_table,
    unordered_map<string, net_map> &net_table,
    unordered_map<string, vector<string>> &inst2nets,
    unordered_map<string, pin_map> &pin_table,
    unordered_map<string, int> &width_dbu,
    int k_neighbors,
    long long &total_hpwl,
    unordered_map<string, long long> &net2hpwl) // ★多這個
{
    bool any_improved = false;
    const int ROWS = (int)rowstripes.size();
    if (ROWS <= 1)
        return false;

    for (int r = 0; r < ROWS; ++r)
    {
        auto &row = rowstripes[r];
        const int N = (int)row.placed_insts.size();
        if (N == 0)
            continue;

        for (int ia = 0; ia < N; ++ia)
        {
            inst_info &A = row.placed_insts[ia];
            if (A.mov_or_fix != "PLACED")
                continue; // ★避免動 fixed cell

            // 只在 A 的最佳區域「不在本列」時才嘗試
            rect region = compute_optimal_region_for_inst(
                A.id, net_table, pin_table, inst_table, inst2nets);

            // row 高度的中線
            int row_mid_y = row.y_origin_dbu; // 如果你有 row_height，就 + row_height/2
            int reg_mid_y = (region.y_bottom + region.y_top) / 2;

            // 若 row 中線落在 optimal region 裡，就不用 vertical swap
            if (reg_mid_y >= row.y_origin_dbu &&
                reg_mid_y <= row.y_origin_dbu)
            {
                // （如果有 row 高度要改成 [y_origin, y_origin+height]）
                continue;
            }

            // 依目標方向決定掃描順序：優先往 optimal region 的方向
            vector<int> dr_order;
            if (reg_mid_y > row_mid_y)
                dr_order = {+1, -1, +2, -2};
            else
                dr_order = {-1, +1, -2, +2};

            // 先收候選：用 A.x±W 的寬鬆窗，並以 "y向收益" 粗估排序
            const int W = 40 * row.x_spacing_dbu; // 可照需要調
            const int y_med = (region.y_bottom + region.y_top) / 2;

            vector<tuple<int, int, int>> pool; // (approx_gain, target_row, ib)

            for (int dr : dr_order)
            {
                int tr = r + dr;
                if (tr < 0 || tr >= ROWS)
                    continue;
                auto &rowT = rowstripes[tr];

                for (int ib = 0; ib < (int)rowT.placed_insts.size(); ++ib)
                {
                    const auto &B = rowT.placed_insts[ib];
                    if (B.id == A.id)
                        continue;
                    if (B.mov_or_fix != "PLACED")
                        continue; // ★不要動 fixed

                    // 放寬 x 視窗：只要求 B.x 靠近 A.x
                    if (B.x_dbu < A.x_dbu - W || B.x_dbu > A.x_dbu + W)
                        continue;

                    // 粗估 y 向收益（越大越好）：把 A 從 row.y 移到 rowT.y
                    int approx_gain = std::abs(A.y_dbu - y_med) - std::abs(rowT.y_origin_dbu - y_med);

                    pool.emplace_back(approx_gain, tr, ib);
                }
            }

            if (pool.empty())
                continue;

            // 取前 k_neighbors 個最高的 approx_gain
            if ((int)pool.size() > k_neighbors)
            {
                std::nth_element(pool.begin(),
                                 pool.begin() + k_neighbors,
                                 pool.end(),
                                 [](const auto &p, const auto &q)
                                 { return std::get<0>(p) > std::get<0>(q); });
                pool.resize(k_neighbors);
            }
            else
            {
                std::sort(pool.begin(), pool.end(),
                          [](const auto &p, const auto &q)
                          { return std::get<0>(p) > std::get<0>(q); });
            }

            long long best_delta = 0;
            int best_tr = -1, best_ib = -1;

            // 逐一檢 legality + 真正 ΔHPWL
            for (auto &cand : pool)
            {
                int tr = std::get<1>(cand);
                int ib = std::get<2>(cand);
                auto &rowT = rowstripes[tr];

                // ★用 fast 版 legality，會檢 row_blocks
                bool legal = check_legal_swap_two_rows_fast(
                    row, ia, rowT, ib, width_dbu);
                if (!legal)
                    continue;

                const auto &B = rowT.placed_insts[ib];

                int ax = B.x_dbu;
                int ay = rowT.y_origin_dbu;
                int bx = A.x_dbu;
                int by = row.y_origin_dbu;

                // ★改成 delta_hpwl_cached，用 net2hpwl
                long long d = delta_hpwl_cached(
                    A.id, B.id,
                    inst_table, net_table, inst2nets, pin_table,
                    net2hpwl, ax, ay, bx, by);

                if (d < best_delta)
                {
                    best_delta = d;
                    best_tr = tr;
                    best_ib = ib;
                }
            }

            if (best_ib >= 0 && best_delta < 0)
            {
                any_improved = true;
                total_hpwl += best_delta;

                auto &rowT = rowstripes[best_tr];
                auto &B = rowT.placed_insts[best_ib];

                // 更新 inst_table
                int Ax = B.x_dbu, Ay = B.y_dbu;
                int Bx = A.x_dbu, By = A.y_dbu;

                inst_table[A.id].x_dbu = Ax;
                inst_table[A.id].y_dbu = Ay;
                inst_table[B.id].x_dbu = Bx;
                inst_table[B.id].y_dbu = By;

                string tmp_orient = inst_table[A.id].orient;
                inst_table[A.id].orient = inst_table[B.id].orient;
                inst_table[B.id].orient = tmp_orient;

                inst_info a_new = inst_table[A.id];
                inst_info b_new = inst_table[B.id];

                row.placed_insts[ia] = b_new;
                rowT.placed_insts[best_ib] = a_new;

                // ★跟 GS 一樣，重算 A/B nets 的 HPWL 更新 net2hpwl（時間很小）
                for (const auto &netname : inst2nets[A.id])
                {
                    long long newv = calculate_HPWL(
                        net_table[netname].connected_insts,
                        inst_table, pin_table);
                    net2hpwl[netname] = newv;
                }
                for (const auto &netname : inst2nets[B.id])
                {
                    long long newv = calculate_HPWL(
                        net_table[netname].connected_insts,
                        inst_table, pin_table);
                    net2hpwl[netname] = newv;
                }

                // ★保持你原本註解：row 之後再整體 sort-by-x，一次做
                // （建議在整個 vertical_swap phase 結束後，對所有 row 做一次
                //   sort(row.placed_insts.begin(), row.placed_insts.end(),
                //        [](const inst_info&a,const inst_info&b){return a.x_dbu<b.x_dbu;});）
            }
        }
    }
    return any_improved;
}

// 假設 row 已按 x 排好
static bool local_reorder_w3(
    RowStripe &row,
    unordered_map<string, int> &width_dbu,
    unordered_map<string, inst_info> &inst_table,
    unordered_map<string, net_map> &net_table,
    unordered_map<string, vector<string>> &inst2nets,
    unordered_map<string, pin_map> &pin_table,
    long long &total_hpwl)
{
    bool improved = false;
    if (row.placed_insts.size() < 3)
        return false;

    const int row_xL = row.x_origin_dbu;
    const int row_xR = row.x_origin_dbu + row.x_sites * row.x_spacing_dbu;

    // 所有 3! 排列
    const int perms[6][3] = {
        {0, 1, 2},
        {0, 2, 1},
        {1, 0, 2},
        {1, 2, 0},
        {2, 0, 1},
        {2, 1, 0}};

    for (int i = 0; i + 2 < (int)row.placed_insts.size(); ++i)
    {
        // 取出 window 的三顆 cell
        inst_info c0 = row.placed_insts[i];
        inst_info c1 = row.placed_insts[i + 1];
        inst_info c2 = row.placed_insts[i + 2];

        // 若有任何 type 找不到寬度就跳過
        if (!width_dbu.count(c0.type) ||
            !width_dbu.count(c1.type) ||
            !width_dbu.count(c2.type))
            continue;

        int w0 = width_dbu[c0.type];
        int w1 = width_dbu[c1.type];
        int w2 = width_dbu[c2.type];

        // window 左邊界 (固定) ＋ 目前三顆的最右邊
        int x_start = c0.x_dbu;
        int old_end = max({c0.x_dbu + w0,
                           c1.x_dbu + w1,
                           c2.x_dbu + w2});

        // 收集受影響的 net
        unordered_set<string> affected_nets;
        for (auto &id : {c0.id, c1.id, c2.id})
        {
            auto it = inst2nets.find(id);
            if (it != inst2nets.end())
            {
                affected_nets.insert(it->second.begin(), it->second.end());
            }
        }

        if (affected_nets.empty())
            continue;

        // 計算原本的 cost
        long long old_cost = 0;
        for (const auto &netname : affected_nets)
        {
            auto &net = net_table[netname];
            old_cost += calculate_HPWL(net.connected_insts, inst_table, pin_table);
        }

        long long best_delta = 0;
        int best_perm_idx = -1;

        // 記錄原始 x，以便試完再還原
        int orig_x0 = inst_table[c0.id].x_dbu;
        int orig_x1 = inst_table[c1.id].x_dbu;
        int orig_x2 = inst_table[c2.id].x_dbu;

        // 嘗試六種排列
        for (int p = 0; p < 6; ++p)
        {
            int order0 = perms[p][0];
            int order1 = perms[p][1];
            int order2 = perms[p][2];

            // 從 x_start 往右「緊密排」，但不改變 window 左邊界
            int x_cur = x_start;

            auto place_one = [&](inst_info &c, int w) -> pair<int, int>
            {
                int L = x_cur;
                int R = L + w;
                x_cur = R;
                return {L, R};
            };

            inst_info *cells[3] = {&c0, &c1, &c2};
            int widths[3] = {w0, w1, w2};

            // 依序算新位置
            auto [L0, R0] = place_one(*cells[order0], widths[order0]);
            auto [L1, R1] = place_one(*cells[order1], widths[order1]);
            auto [L2, R2] = place_one(*cells[order2], widths[order2]);

            // 不能超出 row 邊界
            if (!inside_row_usable_range(L0, R0, row_xL, row_xR) ||
                !inside_row_usable_range(L1, R1, row_xL, row_xR) ||
                !inside_row_usable_range(L2, R2, row_xL, row_xR))
                continue;

            // 檢查不壓到 blocks
            if (!no_overlap_with_blocks(L0, R0, row.row_blocks) ||
                !no_overlap_with_blocks(L1, R1, row.row_blocks) ||
                !no_overlap_with_blocks(L2, R2, row.row_blocks))
                continue;

            // 注意：我們把 window 壓緊到左邊，右界 x_cur <= old_end
            // 因此不會撞到右邊的第 i+3 顆 cell

            // 暫時更新 inst_table 中的 x (不動 y)
            inst_table[cells[order0]->id].x_dbu = L0;
            inst_table[cells[order1]->id].x_dbu = L1;
            inst_table[cells[order2]->id].x_dbu = L2;

            // 計算新 cost
            long long new_cost = 0;
            for (const auto &netname : affected_nets)
            {
                auto &net = net_table[netname];
                new_cost += calculate_HPWL(net.connected_insts, inst_table, pin_table);
            }

            long long delta = new_cost - old_cost;
            if (delta < best_delta)
            {
                best_delta = delta;
                best_perm_idx = p;
            }

            // 還原 inst_table 中的 x
            inst_table[c0.id].x_dbu = orig_x0;
            inst_table[c1.id].x_dbu = orig_x1;
            inst_table[c2.id].x_dbu = orig_x2;
        }

        if (best_perm_idx >= 0 && best_delta < 0)
        {
            // 真的採用 best 排列
            int order0 = perms[best_perm_idx][0];
            int order1 = perms[best_perm_idx][1];
            int order2 = perms[best_perm_idx][2];

            inst_info *cells[3] = {&c0, &c1, &c2};
            int widths[3] = {w0, w1, w2};

            int x_cur = x_start;
            auto place_one = [&](inst_info &c, int w) -> int
            {
                int L = x_cur;
                x_cur = L + w;
                return L;
            };

            int new_x0 = place_one(*cells[order0], widths[order0]);
            int new_x1 = place_one(*cells[order1], widths[order1]);
            int new_x2 = place_one(*cells[order2], widths[order2]);

            // 寫回 inst_table
            inst_table[cells[order0]->id].x_dbu = new_x0;
            inst_table[cells[order1]->id].x_dbu = new_x1;
            inst_table[cells[order2]->id].x_dbu = new_x2;

            // row 的順序就按照 x 排序就好
            row.placed_insts[i] = inst_table[cells[order0]->id];
            row.placed_insts[i + 1] = inst_table[cells[order1]->id];
            row.placed_insts[i + 2] = inst_table[cells[order2]->id];

            // 保險：整列再 sort 一次
            std::sort(row.placed_insts.begin(), row.placed_insts.end(),
                      [](const inst_info &a, const inst_info &b)
                      {
                          return a.x_dbu < b.x_dbu;
                      });

            total_hpwl += best_delta;
            improved = true;
        }
    }

    return improved;
}

//---------------- Parsing ----------------//
void parseLEF(ifstream &lef, unordered_map<string, macro_info> &library_macros)
{
    cout << "Reading LEF...\n";
    string line;
    while (getline(lef, line))
    {
        std::istringstream lef_info(line);
        string word;
        lef_info >> word;

        if (word == "DATABASE")
        {
            string dbu_per_micron;
            lef_info >> word;
            lef_info >> dbu_per_micron;
            dbu = stoi(dbu_per_micron);
            cout << "DBU per micron: " << dbu_per_micron << endl;
        }

        if (word == "SITE")
        {
            while (getline(lef, line))
            {
                if (line.find("SIZE") != string::npos)
                {
                    std::istringstream size_info(line);
                    string skip;
                    string X_size_micron, Y_size_micron;
                    size_info >> skip >> X_size_micron >> skip >> Y_size_micron;
                    cout << "Site Size: " << X_size_micron << " x " << Y_size_micron << endl;
                    break;
                }
            }
        }
        if (word == "MACRO")
        {
            string macro_name;
            lef_info >> macro_name;
            macro_info current_macro;
            while (getline(lef, line))
            {
                if (line.find("CLASS") != string::npos)
                {
                    std::istringstream class_info(line);
                    string skip, macro_class;
                    class_info >> skip >> macro_class;
                    current_macro.class_type = macro_class;
                }
                else if (line.find("SIZE") != string::npos)
                {
                    std::istringstream size_info(line);
                    string skip;
                    double X_size, Y_size;
                    size_info >> skip >> X_size >> skip >> Y_size;
                    current_macro.size_x_micron = X_size;
                    current_macro.size_y_micron = Y_size;
                }
                else if (line.find("SITE") != string::npos)
                {
                    std::istringstream site(line);
                    string skip, site_name;
                    site >> skip >> site_name;
                    current_macro.site = site_name;
                }
                if (line.find("END") != string::npos && line.find(macro_name) != string::npos)
                {
                    break;
                }
            }
            library_macros.insert({macro_name, current_macro});
        }
    }
}

void parseDEF(ifstream &def, vector<row_info> &row_table,
              unordered_map<string, inst_info> &inst_table,
              unordered_map<string, pin_map> &pin_table,
              unordered_map<string, net_map> &net_table,
              unordered_map<string, inst_info> &block_map)
{
    cout << "Reading DEF...\n";
    string line;
    while (getline(def, line))
    {
        std::istringstream def_info(line);
        string word;
        def_info >> word;
        if (word == "UNITS")
        {
            string dbu_per_micron;
            def_info >> word >> word >> dbu_per_micron;
            cout << "DBU per micron: " << dbu_per_micron << endl;
        }
        if (word == "ROW")
        {
            string row_name, site, origin_x_dbu, origin_y_dbu, orient, x_site_size, y_site_size, x_spacing_dbu, y_spacing_dbu, skip;
            def_info >> row_name >> site >> origin_x_dbu >> origin_y_dbu >> orient >> skip >> x_site_size >> skip >> y_site_size >> skip >> x_spacing_dbu >> y_spacing_dbu;
            row_info current_row;
            current_row.row_name = row_name;
            current_row.site = site;
            current_row.origin_x_dbu = stoi(origin_x_dbu);
            current_row.origin_y_dbu = stoi(origin_y_dbu);
            current_row.orient = orient;
            current_row.x_site_size = stoi(x_site_size);
            current_row.y_site_size = stoi(y_site_size);
            current_row.x_spacing_dbu = stoi(x_spacing_dbu);
            // site_width_dbu = current_row.x_spacing_dbu;
            row_table.push_back(current_row);
        }
        if (word == "COMPONENTS")
        {
            int num_components;
            def_info >> num_components;
            for (int i = num_components; i >= 1; i--)
            {
                getline(def, line);
                std::istringstream comp_info(line);
                string skip, comp_num, comp_type;
                inst_info current_inst;
                comp_info >> skip >> comp_num >> comp_type;
                current_inst.id = comp_num;
                current_inst.type = comp_type;
                int count = 0;
                while (line.find("+") == string::npos)
                {
                    getline(def, line);
                    count++;
                }
                string mov_or_fix, x_dbu, y_dbu, orient;
                if (count == 0)
                {
                    comp_info >> skip >> mov_or_fix >> skip >> x_dbu >> y_dbu >> skip >> orient;
                }
                else
                {
                    std::istringstream place_info(line);
                    place_info >> skip >> mov_or_fix >> skip >> x_dbu >> y_dbu >> skip >> orient;
                }
                while (line.find(";") == string::npos)
                {
                    getline(def, line);
                }
                current_inst.mov_or_fix = mov_or_fix;
                current_inst.x_dbu = stoi(x_dbu);
                current_inst.y_dbu = stoi(y_dbu);
                current_inst.orient = orient;
                if (mov_or_fix == "FIXED")
                {
                    block_map.insert({comp_num, current_inst});
                }

                inst_table.insert({comp_num, current_inst});
            }
        }
        if (word == "PINS")
        {
            int num_pins;
            def_info >> num_pins;
            for (int i = 0; i < num_pins; i++)
            {
                getline(def, line);
                std::istringstream pin_info(line);
                pin_map current_pin;
                string pin_name, pin_in_net, pin_x_dbu, pin_y_dbu, orient, skip;
                pin_info >> skip >> pin_name >> pin_in_net;
                current_pin.net_name = pin_in_net;
                while (line.find("PLACED") == string::npos)
                {
                    getline(def, line);
                }
                std::istringstream placed_info(line);
                placed_info >> skip >> skip >> skip >> pin_x_dbu >> pin_y_dbu >> skip >> orient;
                current_pin.x_dbu = stoi(pin_x_dbu);
                current_pin.y_dbu = stoi(pin_y_dbu);
                current_pin.orient = orient;
                pin_table.insert({pin_name, current_pin});
            }
        }
        if (word == "NETS")
        {
            int num_nets;
            def_info >> num_nets;
            for (int i = 0; i < num_nets; i++)
            {
                getline(def, line);
                std::istringstream net_info(line);
                string skip, net_name;
                net_map current_net;
                current_net.connected_insts.clear();
                net_info >> skip >> net_name;

                getline(def, line);
                string tmp;
                while (line.find(";") == string::npos)
                {
                    tmp += " " + line;
                    getline(def, line);
                }
                tmp += ";";
                for (char &c : tmp)
                {
                    if (c == '(' || c == ')')
                        c = ' ';
                }
                std::istringstream connect_info(tmp);
                string insts, part;
                do
                {
                    connect_info >> insts >> part;
                    if (insts == "PIN")
                        insts = part;
                    if (insts == ";")
                        break;
                    current_net.connected_insts.push_back(insts);
                } while (true);
                net_table.insert({net_name, current_net});
            }
        }
    }
}

//-------------------------verify correctness helper-------------------------------------//

static inline void push_err(std::vector<std::string> *all, const std::string &msg)
{
    if (all)
        all->push_back(msg);
}

bool verify_rows_verbose(
    const std::vector<RowStripe> &rows,
    const std::unordered_map<std::string, macro_info> &lib,
    std::string *first_error = nullptr,
    std::vector<std::string> *all_errors = nullptr)
{
    auto fail_once = [&](const std::string &m)
    {
        if (first_error && first_error->empty())
            *first_error = m;
        push_err(all_errors, m);
        return false;
    };

    if (dbu <= 0)
    {
        return fail_once("[GLOBAL] dbu <= 0");
    }

    bool ok = true;

    for (size_t r = 0; r < rows.size(); ++r)
    {

        const auto &row = rows[r];
        if (row.x_spacing_dbu <= 0)
        {
            return fail_once("[GLOBAL] site_width_dbu <= 0");
        }
        const long long row_xL = static_cast<long long>(row.x_origin_dbu);
        const long long row_xR = row_xL + 1LL * row.x_sites * row.x_spacing_dbu;

        long long prev_end = row_xL; // 以 x 升冪假設檢查「與前一顆相交」

        // 基本 row 邏輯檢查
        if (row_xR <= row_xL)
        {
            ok = false;
            std::ostringstream ss;
            ss << "[ROW] r=" << r << " invalid row bound: [" << row_xL << "," << row_xR << ")";
            fail_once(ss.str());
            continue; // 這列都無效，跳過細項
        }

        for (size_t i = 0; i < row.placed_insts.size(); ++i)
        {
            const auto &inst = row.placed_insts[i];

            // 取 macro width
            auto it = lib.find(inst.type);
            if (it == lib.end())
            {
                ok = false;
                std::ostringstream ss;
                ss << "[LIB] r=" << r << " i=" << i << " id=" << inst.id
                   << " unknown macro type: " << inst.type;
                fail_once(ss.str());
                continue;
            }
            const long long width_dbu = static_cast<long long>(std::llround(it->second.size_x_micron * dbu));

            // site 對齊
            if (inst.x_dbu % row.x_spacing_dbu != 0)
            {
                ok = false;
                std::ostringstream ss;
                ss << "[ALIGN] r=" << r << " i=" << i << " id=" << inst.id
                   << " x=" << inst.x_dbu << " not aligned to site_width_dbu=" << row.x_spacing_dbu;
                fail_once(ss.str());
            }

            const long long xL = static_cast<long long>(inst.x_dbu);
            const long long xR = xL + width_dbu;

            // row 邊界
            if (xL < row_xL || xR > row_xR)
            {
                ok = false;
                std::ostringstream ss;
                ss << "[BOUND] r=" << r << " i=" << i << " id=" << inst.id
                   << " seg=[" << xL << "," << xR << ") out of row=[" << row_xL << "," << row_xR << ")";
                fail_once(ss.str());
            }

            // 與前一顆不重疊（允許相切）
            if (i > 0 && xL < prev_end)
            {
                ok = false;
                std::ostringstream ss;
                ss << "[OVERLAP-PREV] r=" << r << " i=" << i << " id=" << inst.id
                   << " xL=" << xL << " < prev_end=" << prev_end;
                fail_once(ss.str());
            }

            // 與下一顆不重疊（允許相切）
            if (i + 1 < row.placed_insts.size())
            {
                const long long next_x = static_cast<long long>(row.placed_insts[i + 1].x_dbu);
                if (xR > next_x)
                {
                    ok = false;
                    std::ostringstream ss;
                    ss << "[OVERLAP-NEXT] r=" << r << " i=" << i << " id=" << inst.id
                       << " xR=" << xR << " > next_x=" << next_x << " (next id=" << row.placed_insts[i + 1].id << ")";
                    fail_once(ss.str());
                }
            }

            // 與 blocks 不相交（假設 block 為半開區間 [bL,bR)）
            for (size_t bi = 0; bi < row.row_blocks.size(); ++bi)
            {
                const long long bL = static_cast<long long>(row.row_blocks[bi].first);
                const long long bR = static_cast<long long>(row.row_blocks[bi].second);
                if (xL < bR && xR > bL)
                {
                    ok = false;
                    std::ostringstream ss;
                    ss << "[BLOCK] r=" << r << " i=" << i << " id=" << inst.id
                       << " seg=[" << xL << "," << xR << ") overlaps block[" << bi << "]="
                       << "[" << bL << "," << bR << ")";
                    fail_once(ss.str());
                }
            }

            prev_end = std::max(prev_end, xR);
        }
    }

    return ok && (first_error ? first_error->empty() : true);
}

//------------------------write output-------------------------------------//
// 去掉左邊空白的小工具
static std::string ltrim(const std::string &s)
{
    size_t i = 0;
    while (i < s.size() && std::isspace((unsigned char)s[i]))
        ++i;
    return s.substr(i);
}

#include <filesystem>
#include <cerrno>
#include <cstring>

namespace fs = std::filesystem;

void rewrite_def_coords_only(
    std::string &in_def,
    std::string &out_def,
    std::unordered_map<std::string, inst_info> &inst_table,
    std::vector<row_info> &row_table)
{
    // 1) 確保 output 資料夾存在
    try
    {
        fs::path out_path(out_def);
        fs::path parent = out_path.parent_path();
        if (!parent.empty() && !fs::exists(parent))
        {
            fs::create_directories(parent);
            std::cout << "[INFO] Created output dir: " << parent << "\n";
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "[ERROR] create_directories failed for \"" << out_def
                  << "\": " << e.what() << "\n";
        // 先不要 return，萬一目錄本來就存在
    }

    // 2) 開 input DEF
    std::ifstream in(in_def);
    if (!in.is_open())
    {
        std::cerr << "[ERROR] Cannot open input DEF: " << in_def
                  << " (errno=" << errno << " " << std::strerror(errno) << ")\n";
        return;
    }

    // 3) 開 output DEF（不存在就會新建，存在就覆寫）
    std::ofstream out(out_def); // 預設是 trunc 模式
    if (!out.is_open())
    {
        std::cerr << "[ERROR] Cannot open output DEF: " << out_def
                  << " (errno=" << errno << " " << std::strerror(errno) << ")\n";
        return;
    }

    std::string line;
    bool in_components = false;
    std::string current_inst_name;
    long long updated_cnt = 0;

    while (std::getline(in, line))
    {
        std::string trimmed = ltrim(line);

        // 進入 / 離開 COMPONENTS
        if (trimmed.rfind("COMPONENTS", 0) == 0)
        {
            in_components = true;
            current_inst_name.clear();
            out << line << "\n";
            continue;
        }
        if (trimmed.rfind("END COMPONENTS", 0) == 0)
        {
            in_components = false;
            current_inst_name.clear();
            out << line << "\n";
            continue;
        }

        if (in_components)
        {
            // 例: "  - U1234 NAND2_X1"
            if (trimmed.rfind("- ", 0) == 0)
            {
                std::istringstream iss(trimmed);
                std::string dash, inst_name;
                iss >> dash >> inst_name;
                current_inst_name = inst_name;

                out << line << "\n"; // header 行照抄
                continue;
            }

            // 例: "    + PLACED ( 100 200 ) N ;"
            if (!current_inst_name.empty() && trimmed.rfind("+", 0) == 0)
            {
                auto it = inst_table.find(current_inst_name);
                if (it != inst_table.end())
                {
                    const inst_info &inst = it->second;

                    out << "    + " << inst.mov_or_fix
                        << " ( " << inst.x_dbu << " " << inst.y_dbu << " ) "
                        << inst.orient << " ;\n";

                    ++updated_cnt;
                    continue; // 不寫原本那行
                }
            }
        }

        // 其他行全部照抄
        out << line << "\n";
    }

    in.close();
    out.close();

    std::cout << "[INFO] rewrite_def_coords_only: updated "
              << updated_cnt << " components for " << out_def << "\n";
}

//---------------------------------main-----------------------//

int main(int argc, char *argv[])
{
    using namespace std::chrono;
    //----------計時開始----------//
    auto start = high_resolution_clock::now(); // 開始計時
    //--------------------------//
    // Expect lef and def
    // freopen("log.txt", "w", stdout);
    if (argc < 4)
    {
        cerr << "Usage: " << argv[0] << " <lef_file> <def_file>" << endl;
        return 1;
    }

    ifstream lef(argv[1]);
    if (!lef.is_open())
    {
        cerr << "Error: Cannot open file " << argv[1] << endl;
        return 1;
    }
    ifstream def(argv[2]);
    if (!def.is_open())
    {
        cerr << "Error: Cannot open file " << argv[2] << endl;
        return 1;
    }
    ifstream def_output(argv[3]);

    unordered_map<string, macro_info> library_macros;
    vector<row_info> row_table;
    unordered_map<int, int> y2row_idx;
    unordered_map<string, inst_info> inst_table;
    unordered_map<string, pin_map> pin_table;
    unordered_map<string, net_map> net_table;
    unordered_map<string, inst_info> block_map;

    parseLEF(lef, library_macros);
    parseDEF(def, row_table, inst_table, pin_table, net_table, block_map);
    auto now = chrono::high_resolution_clock::now();
    cout << " Time: " << chrono::duration_cast<chrono::milliseconds>(now - start).count() << endl;

    vector<RowStripe> row_stripes = make_cell_map(row_table, inst_table, block_map);
    calculate_block_bound(block_map, library_macros, row_stripes);

    unordered_map<string, vector<string>> inst2nets = build_inst2nets(net_table);
    unordered_map<string, int> width_dbu = build_macro_width_dbu(library_macros);
    makey2row_idx(row_stripes, y2row_idx);
    now = chrono::high_resolution_clock::now();

    // Diagnostics: missing width types
    /*size_t miss = 0;
    for (auto &kv : inst_table)
    {
        if (!width_dbu.count(kv.second.type))
            miss++;
    }
    std::string first_err;
    std::vector<std::string> all_errs;

    /*bool ok = verify_rows_verbose(row_stripes, library_macros, &first_err, &all_errs);

    if (!ok) {
        // 印第一個錯誤
        cout<<"error first"<<endl;
        std::cerr << "Verify failed: " << first_err << "\n";
        // 或印全部錯誤
        for (auto& s : all_errs) std::cerr << s << "\n";
    }
    else{
        cout<<"OK first"<<endl;
    }*/

    // cout<<"Initial table"<<endl;
    auto net2hpwl = build_net2hpwl_cache(net_table, inst_table, pin_table);
    long long total_hpwl = 0;
    for (auto &kv : net2hpwl)
        total_hpwl += kv.second;

    cout << "Initial HPWL: " << total_hpwl << endl;
    int x_w = 0;
    int max_row_delta = 0;
    int K_near_for_GS = 0;
    int vs_k_neighbors = 0;    // 原本 20，建議 50~80 視效能
    int max_outer_iters = 0;   // 外圈最多 12~15 輪
    int no_improve_rounds = 0; // 連續無改善的輪數，做早停
    bool flag;
    bool flag1;
    // 調參
    if (inst_table.size() > 200000 && inst_table.size() < 600000)
    {
        x_w = 20;
        K_near_for_GS = 8;
        max_row_delta = 2;
        vs_k_neighbors = 10;
        max_outer_iters = 1;
        TIME_LIMIT_MS = 200000;
        flag = true;
        flag1 = true;
    }
    else if (inst_table.size() >= 600000)
    {
        x_w = 20;
        K_near_for_GS = 3;
        max_row_delta = 3;
        vs_k_neighbors = 5;
        max_outer_iters = 1;
        TIME_LIMIT_MS = 230000;
        flag = true;
        flag1 = true;
    }
    else
    {
        x_w = 80;
        K_near_for_GS = 60;
        max_row_delta = 10;
        vs_k_neighbors = 40;
        max_outer_iters = 8;
        TIME_LIMIT_MS = 270000;
        flag = true;
        flag1 = true;
    }

    // print_row_with_cells(row_stripes);
    reset_global_swap_timers();
    for (int iter = 0; iter < max_outer_iters; ++iter)
    {
        bool improved_round = false;
        auto now = chrono::high_resolution_clock::now();
        if (chrono::duration_cast<chrono::milliseconds>(now - start).count() > TIME_LIMIT_MS)
        {
            cout << "Time limit reached, stop optimization.\n";
            cout << " Time: " << chrono::duration_cast<chrono::milliseconds>(now - start).count() << endl;
            break;
        }

        // 1) Global Swap
            bool imp_gs = global_swap_one_pass(
                row_stripes, inst_table, net_table, inst2nets,
                pin_table, width_dbu, max_row_delta, total_hpwl,
                net2hpwl, K_near_for_GS, x_w, y2row_idx);
            improved_round |= imp_gs;
            cout << "[Iter " << iter << "][GS ] HPWL = " << total_hpwl << (imp_gs ? " (improved)" : "") << "\n";
            now = chrono::high_resolution_clock::now();
            if (chrono::duration_cast<chrono::milliseconds>(now - start).count() > TIME_LIMIT_MS)
            {
                cout << "Time limit reached, stop optimization.\n";
                cout << " Time: " << chrono::duration_cast<chrono::milliseconds>(now - start).count() << endl;
                break;
            }

    

            // 2) Vertical Swap
            bool imp_vs = vertical_swap_one_pass(
                row_stripes, inst_table, net_table, inst2nets,
                pin_table, width_dbu, vs_k_neighbors, total_hpwl, net2hpwl);
            improved_round |= imp_vs;
            cout << "[Iter " << iter << "][VS ] HPWL = " << total_hpwl << (imp_vs ? " (improved)" : "") << "\n";
            now = chrono::high_resolution_clock::now();
            if (chrono::duration_cast<chrono::milliseconds>(now - start).count() > TIME_LIMIT_MS)
            {
                cout << "Time limit reached, stop optimization.\n";
                cout << " Time: " << chrono::duration_cast<chrono::milliseconds>(now - start).count() << endl;
                break;
            }
        // 3) Local Reordering（window-3）
        bool imp_lr = false;
        for (auto &row : row_stripes)
        {
            bool row_improved = local_reorder_w3(
                row, width_dbu, inst_table, net_table, inst2nets, pin_table, total_hpwl);
            imp_lr |= row_improved;
        }
        improved_round |= imp_lr;
        cout << "[Iter " << iter << "][LR ] HPWL = " << total_hpwl << (imp_lr ? " (improved)" : "") << "\n";
        now = chrono::high_resolution_clock::now();
        if (chrono::duration_cast<chrono::milliseconds>(now - start).count() > TIME_LIMIT_MS)
        {
            cout << "Time limit reached, stop optimization.\n";
            cout << " Time: " << chrono::duration_cast<chrono::milliseconds>(now - start).count() << endl;
            break;
        }
        net2hpwl = build_net2hpwl_cache(net_table, inst_table, pin_table);
        total_hpwl = 0;
        for (auto &kv : net2hpwl)
            total_hpwl += kv.second;
        // cout<<total_hpwl<<endl;
        //  早停條件：連續兩輪無改善就停
        if (!improved_round)
        {
            ++no_improve_rounds;
            if (no_improve_rounds >= 2)
                break;
        }
        else
        {
            no_improve_rounds = 0;
        }
        now = chrono::high_resolution_clock::now();
        if (chrono::duration_cast<chrono::milliseconds>(now - start).count() > TIME_LIMIT_MS)
        {
            cout << "Time limit reached, stop optimization.\n";
            break;
        }
    }

    print_global_swap_timers();
    // Refresh cache & recompute for consistency
    // cout<<"table 2"<<endl;
    net2hpwl = build_net2hpwl_cache(net_table, inst_table, pin_table);
    /*for(auto &kx :temp){
        for(auto &kv:net2hpwl){
            if(kx.first==kv.first){
                if(kx.second!=kv.second){
                    cout<<kx.first<<" "<<kx.second<<" "<<kv.first<<" "<<kv.second<<endl;
                }
                else{
                    //cout<<"no differ"<<endl;
                }
            }
        }
    }*/
    total_hpwl = 0;
    for (auto &kv : net2hpwl)
        total_hpwl += kv.second;
    // print_row_with_cells(row_stripes);
    cout << "Final HPWL: " << total_hpwl << endl;
    /*auto ok = verify_rows_verbose(row_stripes, library_macros, &first_err, &all_errs);

    if (!ok)
    {
        // 印第一個錯誤
        cout << "error last" << endl;
        std::cerr << "Verify failed: " << first_err << "\n";
        // 或印全部錯誤
        for (auto &s : all_errs)
            std::cerr << s << "\n";
    }
    else
    {
        cout << "OK last" << endl;
    }*/
    string input_def = argv[2];
    string output_def = argv[3];
    rewrite_def_coords_only(input_def, output_def, inst_table, row_table);

    //----------記錄時間----------//
    auto end = high_resolution_clock::now();                  // 結束計時
    auto duration = duration_cast<milliseconds>(end - start); // 計算時間差
    cout << "Execution time: " << duration.count() << " milliseconds" << endl;
    //--------------------------//
    return 0;
}
