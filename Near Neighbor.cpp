#include <bits/stdc++.h>

using namespace std;

const float EPS = 1e-6;
const float INF = 1e7;

typedef vector<float> Point;

float radius, lambda, RADIUS;
float rho, threshold;
int d, nodes;
vector<Point> points;

float ELLINF(const Point& x, const Point& y) {
  float ans = 0;
  for (int i = 0; i < d; ++i) {
    ans = max(ans, fabs(x[i] - y[i]));
  }
  return ans;
}

struct Node {
  Node *l;
  Node *r;
  bool isSep, isLeaf;
  pair<int, float> hyperplane;
  Point rep_point;
  int ball_point;
  vector<int> leafNodes;
  Node() {
    l = nullptr;
    r = nullptr;
    isSep = isLeaf = false;
  }
};

float L(float beta, float gamma) {
  float res = (log(gamma) / log(beta + gamma)) - 1.0;
  return res;
}

float hyperplane_helper(int coordinate, vector<int>& P) {
  int n = P.size(), j = 0;
  threshold = 1.0 / (4.0 * d);
  sort(P.begin(), P.end(), [&](int x, int y) {
    return points[x][coordinate] < points[y][coordinate];
  });
  float alpha, beta, gamma;
  for (int i = n / 4; i < n; ++i) {
    j = max(i + 1, j);
    while (j < n && points[P[j]][coordinate] < points[P[i]][coordinate] + 2.0 * lambda) ++j;
    alpha = float(i + 1) / n;
    gamma = float(n - j) / n;
    beta = float(j - i - 1) / n;
    if (alpha < threshold || gamma < threshold || log(gamma) > (rho + 1) * log(beta + gamma) + EPS) continue;
    return (points[P[i]][coordinate] + points[P[j]][coordinate]) / 2.0;
  }
  return INF + INF;
}

void cleanse(int n, vector<int>& idx) {
  float perturbation = 1e-5;
  for (int i = 0; i < d; ++i) {
    sort(idx.begin(), idx.end(), [i](int x, int y) {
      return points[x][i] < points[y][i];
    });
    for (int j = 1; j < points.size(); ++j) {
      points[j][i] = max(points[j - 1][i] + perturbation, points[j - 1][i]);
    }
  }
  cerr << "Data Cleaned..\n";
}

pair<int, float> getHyperplane(vector<int>& P) {
  int n = P.size();
  for (int i = 0; i < d; ++i) {
    float dividerPlane = hyperplane_helper(i, P);
    if (dividerPlane > INF) continue;
    return make_pair(i, dividerPlane);
  }
  return make_pair(-1, 0);
}

float getMedian(vector<int>& P, int coordinate) {
  sort(P.begin(), P.end(), [&](int x, int y) {
    return points[x][coordinate] < points[y][coordinate];
  });
  int n = P.size();
  int x = n / 2, y = (n - 1) / 2;
  float res = (points[P[x]][coordinate] + points[P[y]][coordinate]) / 2.0;
  return res;
}

int getBall(vector<int>& P, Point& med, float rad = RADIUS) {
  med.resize(d);
  int C = -1;
  for (int i = 0; i < d; ++i) {
    med[i] = getMedian(P, i);
  }
  vector<int> outball;
  for (int id : P) {
    if (ELLINF(points[id], med) > rad + EPS) {
      outball.push_back(id);
    } else {
      C = id;
    }
  }
  assert(P.size() > outball.size());
  P = outball;
  return C;
}

Node* preprocessTree(vector<int>& P) {
  Node* root = new Node();
  ++nodes;
  if (P.size() <= 20) {
    root->isLeaf = true;
    root->leafNodes = P;
    return root;
  }
  auto divider = getHyperplane(P);
  if (divider.first == -1) {
    root->ball_point = getBall(P, root->rep_point, RADIUS);
    if (!P.empty()) root->l = preprocessTree(P);
    return root;
  }
  int n = P.size();
  vector<int> L, R;
  for (int idx : P) {
    if (points[idx][divider.first] + EPS < divider.second + lambda) {
      L.push_back(idx);
    }
    if (points[idx][divider.first] > divider.second - lambda + EPS) {
      R.push_back(idx);
    }
  }
  if (!(L.size() < n && R.size() < n)) {
    cerr << divider.first << " " << divider.second << "\n";
    assert(0);
  }
  P.clear();
  root->isSep = true;
  root->hyperplane = divider;
  root->l = preprocessTree(L);
  root->r = preprocessTree(R);
  return root;
}

int searchnns(Node* &root, const Point& query) {
  if (root == nullptr) return -1;
  if (root->isLeaf) {
    float dis = INF;
    int res = -1;
    for (int id : root->leafNodes) {
      float cur = ELLINF(query, points[id]);
      if (cur < dis) {
        dis = cur;
        res = id;
      }
    }
    return res;
  }

  if (root->isSep) {
    pair<int, float> hyperplaneData = root->hyperplane;
    if (query[hyperplaneData.first] < hyperplaneData.second) {
      return searchnns(root->l, query);
    }
    else {
      return searchnns(root->r, query);
    }
  }
  if (ELLINF(query, root->rep_point) < RADIUS + lambda + EPS || root->l == nullptr) {
    return root->ball_point;
  }
  else {
    return searchnns(root->l, query);
  }
}

void readData(int n) {
  points.resize(n);
  for (int i = 0; i < n; ++i) {
    points[i].resize(d);
    for (int j = 0; j < d; ++j) {
      cin >> points[i][j];
    }
  }

  cerr << "Read data.\n\n";
}

void readQuery(int q, vector<Point>& QQ) {
  for (int i = 0; i < q; ++i) {
    QQ[i].resize(d);
    for (int j = 0; j < d; ++j) {
      cin >> QQ[i][j];
    }
  }
  cerr << "Queries taken.\n";
}

int main() {
  int n, q;

  ios_base::sync_with_stdio(0);
  cin.tie(0);

  cin >> rho >> threshold >> lambda;
  cin >> n >> q >> d;

  readData(n);
  vector<Point> QQ(q);
  readQuery(q, QQ);

  radius = 2 * ceil(log (log2(4 * d)) / log(1.0 + rho));
  RADIUS = radius * lambda;

  vector<int> idxs(n);
  for (int i = 0; i < n; ++i) idxs[i] = i;
  cleanse(n, idxs);

  Node* root = preprocessTree(idxs);
  cerr << "Done building tree.\n";
  cerr << "Nodes: " << nodes << "\n";

  cout << nodes << "\n" << QQ.size() << "\n" << QQ[0].size() << "\n";

  clock_t start;

  start = clock();
  for (auto pt : QQ) {
    int approx = searchnns(root, pt);
    cout << ELLINF(points[approx], pt) << "\n";
  }
  cout << (clock() - start) / (float) CLOCKS_PER_SEC << "\n";
  return 0;
}
