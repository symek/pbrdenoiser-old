/* 
  Code accompaning:

  A Fresh Look at Generalized Sampling
  Diego Nehab
  Instituto Nacional de Matemática Pura e Aplicada (IMPA)
  Hugues Hoppe
  Microsoft Research 

  http://hhoppe.com/proj/filtering/

*/


#include <iostream>   // std::ostream, std::cout
#include <cmath>      // std::abs(), std::ceil(), std::floor()
#include <algorithm>  // std::min()
#include <vector>     // std::vector<>
#include <array>      // std::array<>
#include <string>     // std::string
#include <cassert>    // assert()
#include <chrono>     // std::chrono for timing benchmarks
using namespace std;
namespace GS {
// Apply to sequence f the inverse discrete convolution given by
// a pre-factored LU decomposition\label{lst:inverse-convolution-begin}
template<size_t M>
void linear_solve(const array<float,M>& L, vector<float>& f) {
  const int m = M, n = int(f.size());
  const float p_inv = 1.f; // Optimized for prescaled kernel where $p=1$
  // Pre-factored decomposition only works when n>m.  Grow sequence f if needed.
  while (f.size() <= m) {
    f.reserve(2*f.size()); // Prevent reallocation during insertions
    f.insert(f.end(), f.rbegin(), f.rend()); // Append a reflection
  }
  int nn = int(f.size()); // New size
  const float L_inf = L[m-1], v_inv = L_inf/(1.f+L_inf);
  // Forward pass: solve $L c'= f$ in-place
  for (int i=1; i<m;  i++)      f[i] -= L[i-1]*f[i-1];
  for (int i=m; i<nn; i++)      f[i] -= L_inf *f[i-1];
  // Reverse pass: solve $U c = c'$ in-place
  f[nn-1] *= p_inv*v_inv;
  for (int i=nn-2; i>=m-1; i--) f[i] = L_inf*(p_inv*f[i]-f[i+1]);
  for (int i=m-2;  i>=0;   i--) f[i] = L[i] *(p_inv*f[i]-f[i+1]);
  f.resize(n); // Truncate back to original size if grown
}/*\label{lst:inverse-convolution-end}*/

// Given sequence f, access f[i] using reflection at boundaries\label{lst:reflect-begin}
static inline float get(const vector<float>& f, int i) {
  return i<0? get(f,-i-1): i>=int(f.size())?
    get(f,2*int(f.size())-i-1): f[i];
}/*\label{lst:reflect-end}*/

// Kernel interface\label{lst:base-kernel-begin}
template<size_t N>
class KernelBase {
  public:
    KernelBase() { b.fill(0.f); }
    // Evaluate kernel at coordinate x
    virtual float operator()(float x) const = 0;
    // Apply the kernel's associated digital filter to sequence f
    virtual void digital_filter(vector<float>& f) const = 0;
    // Kernel support $(\mm\mathtt{support}/2,\mathtt{support}/2]$
    int support() const { return N; }
    // Shift the incremental buffer
    float shift_buffer(float a) {
      rotate(b.begin(), b.begin()+1, b.end());
      swap(b.back(), a);
      return a;
    }
    // Incrementally accumulate a sample into buffer
    virtual void accumulate_buffer(float fu, float u) = 0;
    // Incrementally reconstruct the function from samples in buffer
    virtual float sample_buffer(float u) const = 0;
    virtual string name() const = 0;
    // How much does the kernel integrate to?
    virtual float integral() const { return 1.f; }
  protected:
    array<float,N> b; // Incremental buffer
};/*\label{lst:base-kernel-end}*/

// Simple box kernel\label{lst:box-begin}
struct Box final: KernelBase<1> {
    float operator()(float x) const override {
      return x<=-0.5f || x>0.5f ? 0.f : 1.f;
    }
    void accumulate_buffer(float fu, float) override {
      b[0] += fu;
    }
    float sample_buffer(float) const override { return b[0]; }
    void digital_filter(vector<float>&) const override { }
    string name() const override { return "Box"; }
};/*\label{lst:box-end}*/

// Simple hat kernel\label{lst:hat-begin}
struct Hat final: KernelBase<2> {
    float operator()(float x) const override {
      x = abs(x); return x>1.f ? 0.f : 1.f-x;
    }
    void accumulate_buffer(float fu, float u) override {
      b[0] += fu*(1.f-u); b[1] += fu*u;
    }
    float sample_buffer(float u) const override {
      return b[0]*(1.f-u)+b[1]*u;
    }
    void digital_filter(vector<float>&) const override { }
    string name() const override { return "Hat"; }
};/*\label{lst:hat-end}*/

// Most cubics are $\diff{1}$-continuous, symmetric, and have support 4.
// Factor out common functionality into a class.
template<typename Pieces>
class Symmetric4Pieces: public KernelBase<4> {
  public:
    float operator()(float x) const override final {
      x=abs(x); return x>2.f ? 0.f : x>1.f ? p.k0(2.f-x) : p.k1(1.f-x);
    }
    void accumulate_buffer(float fu, float u) override final {
      b[0] += fu*p.k3(u); b[1] += fu*p.k2(u);
      b[2] += fu*p.k1(u); b[3] += fu*p.k0(u);
    }
    float sample_buffer(float u) const override final {
      return b[0]*p.k3(u)+b[1]*p.k2(u)+b[2]*p.k1(u)+b[3]*p.k0(u);
    }
  private:
    Pieces p; // Polynomial pieces of kernel (k0:[-2,-1], k1:[-1,0], k2:[0,1], k3:[1,2]
};

// Traditional Mitchell-Netravali kernel\label{lst:mitchell-begin}
struct MitchellNetravaliPieces {
  static float k0(float u) {
    return (((7/18.f)*u-1/3.f)*u)*u;
  }
  static float k1(float u) {
    return (((-7/6.f)*u+1.5f)*u+.5f)*u+1/18.f;
  }
  static float k2(float u) {
    return (((7/6.f)*u-2.f)*u)*u+8/9.f;
  }
  static float k3(float u) {
    return (((-7/18.f)*u+5/6.f)*u-.5f)*u+1/18.f;
  }
};

struct MitchellNetravali final:
  Symmetric4Pieces<MitchellNetravaliPieces> {
    void digital_filter(vector<float>&) const override { }
    string name() const override { return "Mitchell-Netravali"; }
};/*\label{lst:mitchell-end}*/


// Traditional Catmull-Rom kernel\label{lst:keys-begin}
struct CatmullRomPieces {
  static float k0(float u) { return ((.5f*u-.5f)*u)*u; }
  static float k1(float u) { return ((-1.5f*u+2.f)*u+.5f)*u; }
  static float k2(float u) { return ((1.5f*u-2.5f)*u)*u+1.f; }
  static float k3(float u) { return ((-.5*u+1.f)*u-.5f)*u; }
};

struct CatmullRom final: Symmetric4Pieces<CatmullRomPieces> {
    void digital_filter(vector<float>&) const override { }
    string name() const override { return "Catmull-Rom"; }
};/*\label{lst:keys-end}*/

// Cubic B-spline kernel pieces (multiplied by 6)
struct Bspline3Pieces {
  static float k0(float u) { return ((u)*u)*u; }
  static float k1(float u) { return ((-3.f*u+3.f)*u+3.f)*u+1.f; }
  static float k2(float u) { return ((3.f*u-6.f)*u)*u+4.f; }
  static float k3(float u) { return ((-u+3.f)*u-3.f)*u+1.f; }
};

// Generalized Cardinal Cubic B-spline kernel\label{lst:cb3-begin}
struct CardinalBspline3 final: Symmetric4Pieces<Bspline3Pieces> {
    void digital_filter(vector<float>& f) const override {
      // Pre-factored L U decomposition of digital filter\label{lst:cb3-prefactored-begin}
      const array<float,8> L{.2f, .26315789f, .26760563f,
        .26792453f, .26794742f, .26794907f, .26794918f, .26794919f};/*\label{lst:cb3-prefactored-end}*/
      linear_solve(L, f);
    }
    string name() const override {
      return "Cardinal Cubic B-spline";
    }
    float integral() const override { return 6.f; }
};/*\label{lst:cb3-end}*/

// Cubic OMOMS kernel pieces (multiplied by 5.25)
struct OMOMS3Pieces {
  static float k0(float u) {
    return ((.875f*u)*u+.125f)*u;
  }
  static float k1(float u) {
    return ((-2.625f*u+2.625f)*u+2.25f)*u+1.f;
  }
  static float k2(float u) {
    return ((2.625f*u-5.25f)*u+.375f)*u+3.25f;
  }
  static float k3(float u) {
    return ((-.875f*u+2.625f)*u-2.75f)*u+1.f;
  }
};

// Generalized Cardinal Cubic O-MOMS3 kernel\label{lst:o3-begin}
struct CardinalOMOMS3 final: Symmetric4Pieces<OMOMS3Pieces> {
    void digital_filter(vector<float>& f) const override {
      // Pre-factored L U decomposition of digital filter\label{lst:o3-prefactored-begin}
      const array<float,9> L{.23529412f, .33170732f, .34266611f,
        .34395774f, .34411062f, .34412872f, .34413087f, .34413112f,
        .34413115f}; /*\label{lst:o3-prefactored-end}*/
      linear_solve(L, f);
    }
    string name() const override { return "Cardinal Cubic OMOMS"; }
    float integral() const override { return 5.25f; }
};/*\label{lst:o3-end}*/

// Simple, intuitive implementation of upsampling\label{lst:intuitive-upsample-begin}
template<typename Kernel> static vector<float>
upsample(vector<float> f, int m, Kernel& k) {
  assert(m >= int(f.size())); // Ensure we are upsampling
  vector<float> g(m); // New sequence of desired size m>f.size()
  k.digital_filter(f); // Apply kernel's associated digital filter
  const float kr = .5f*float(k.support());
  for (int j=0; j<m; j++) { // Index of sample in g
    float x = (j+.5f)/m; // Position in domain [0,1] of both f and g
    float xi = x*f.size()-.5f; // Position in input sequence f
    int il = int(ceil(xi-kr)); // Leftmost sample under kernel support
    int ir = int(floor(xi+kr)); // Rightmost sample under kernel support
    double sum = 0.;
    for (int i=il; i<=ir; i++)
      sum += get(f,i)*k(xi-i);
    g[j] = float(sum);
  }
  return g;
}/*\label{lst:intuitive-upsample-end}*/

// Simple, intuitive implementation of downsampling\label{lst:intuitive-downsample-begin}
template<typename Kernel> static vector<float>
downsample(const vector<float>& f, int m, Kernel& k) {
  assert(m <= int(f.size()));  // Ensure we are downsampling
  float s = float(m)/f.size(); // Scale factor
  const int n = int(f.size());
  const float kr = .5f*float(k.support());
  const bool should_normalize = (f.size()%m != 0);
  vector<float> g(m); // New sequence of desired size m<f.size()
  for (int j=0; j<m; j++) { // Index of sample in g
    float x = (j+.5f)/m; // Position in domain [0,1] of both f and g
    int il = int(ceil ((x-kr/m)*n-.5f)); // Leftmost sample under kernel
    int ir = int(floor((x+kr/m)*n-.5f)); // Rightmost sample under kernel
    if (should_normalize) { // Should normalize?
      double sum = 0., sumw = 0.;  // Sums of values and weights
      for (int i=il; i<=ir; i++) { // Loop over input samples
        float w = k((x-(i+.5f)/n)*m); // Weight for sample
        sum += w*get(f,i); sumw += w; // Accumulate values and weights
      }
      g[j] = k.integral()*float(sum/sumw); // Normalize by summed weights
    } else {
      for (int i=il; i<=ir; i++) {
        g[j] += k((x-(i+.5f)/n)*m)*get(f,i);
      }
      g[j] *= s;
    }
  }
  k.digital_filter(g); // Apply kernel's associated digital filter
  return g;
}/*\label{lst:intuitive-downsample-end}*/

// Advance to next sample
inline bool advance(int& i, int& j, double& u, double inv_s) {
  ++i; u += inv_s;
  if (u < 1.) return false;
  u -= 1.; ++j; return true;
}

// Faster, incremental implementation of upsampling\label{lst:incremental-upsample-begin}
template<typename Kernel> vector<float>
upsample2(vector<float> f, int m, Kernel& k) {
  k.digital_filter(f);
  double inv_s = double(f.size())/m; // Inverse scale factor
  // Output sample position between input samples
  double u = .5*(inv_s+(k.support()+1)%2);
  int fi = -k.support()/2-1;
  assert(f.size() <= m); // Ensure we are upsampling
  vector<float> g(m); // New sequence of desired size m>f.size()
  for (int i=0; i<k.support(); i++) // Initialize incremental buffer
    k.shift_buffer(get(f, ++fi));
  for (int gi=0; gi<m; ) { // Sample reconstruction of f into g[gi]
    g[gi] = k.sample_buffer(u);
    if (advance(gi, fi, u, inv_s))
      k.shift_buffer(get(f,fi));
  }
  return g;
}/*\label{lst:incremental-upsample-end}*/

// Faster, incremental implementation of downsampling\label{lst:incremental-downsample-begin}
template<typename Kernel> vector<float>
downsample2(const vector<float>& f, int m, Kernel& k) {
  const int n = f.size();
  double s = double(m)/n; // Scale factor
  double kr = .5*k.support();
  int fi = int(ceil(((.5-kr)/m)*n-.5));
  // Input sample position between output samples
  double u = ((fi+.5)/n)*m -(.5-kr);
  assert(f.size() >= m); // Ensure we are downsampling
  vector<float> g(m); // New sequence of desired size m>f.size()
  int gi = -k.support();
  for (int i=0; i < k.support(); i++) // Initialize incremental buffer
    k.shift_buffer(0.f);
  while (gi < -1) {
    k.accumulate_buffer(get(f, fi), u);
    if (advance(fi, gi, u, s)) k.shift_buffer(0.f);
  }
  while (1) { // Accumulate weighted f samples into g
    k.accumulate_buffer(get(f, fi), u);
    if (advance(fi, gi, u, s)) {
      if (gi >= m) break;
      g[gi] = s*k.shift_buffer(0.f);
    }
  }
  k.digital_filter(g);
  return g;
}/*\label{lst:incremental-downsample-end}*/

// Output a sequence
static ostream& operator<<(ostream& out, const vector<float>& f) {
  for (int i=0; i<f.size(); i++)
    out << (i+.5)/f.size() << '\t' << f[i] << '\n';
  return out;
}

// Report elapsed lifetime of object
class Timing {
  public:
    Timing(const string& name) : m_name(name), m_start(now()) { }
    ~Timing() {
      using namespace std::chrono;
      double t = duration_cast<microseconds>(now()-m_start).count();
      cout << m_name << " in " << t/1000. << " ms\n";
    }
  private:
    using clock = chrono::high_resolution_clock;
    using time_point = chrono::time_point<clock>;
    static time_point now() { return clock::now(); }
    string m_name;
    time_point m_start;
};

static double
maxerr(const vector<float>& f, const vector<float>& g) {
  double m = 0.;
  assert(f.size() == g.size());
  for (size_t i=0; i<f.size(); i++) {
    m = max(m, abs(double(g[i])-f[i]));
  }
  return m;
}

template<typename Kernel>
void test_performance(int up, int down) {
  Kernel k;
  vector<float> f{0.f, 3.f, 1.f, .5f, 4.f, 2.f}, g(f);
  { Timing t(k.name()+" up");    f = upsample(f, up, k); }
  { Timing t(k.name()+" up2");   g = upsample2(g, up, k); }
  cout << "  " << maxerr(f, g) << "  max error\n";
  { Timing t(k.name()+" down");  f = downsample(f, down, k); }
  { Timing t(k.name()+" down2"); g = downsample2(g, down, k); }
  cout << "  " << maxerr(f, g) << "  max error\n";
}

// Plot kernel impulse response. Pipe through \verb|gnuplot -persist|
template <typename Kernel>
static void
gnuplot(const vector<float>& f, int window, Kernel& k, int w) {
  cout << "set terminal aqua " << window << '\n'; // Change to your needs
  cout << "set title \"" << k.name() << "\"\n";
  cout << "set xrange [" << -.5f*w << ":" << .5f*w << "]\n";
  cout << "plot \"-\" u (" << w << "*($1-0.5)):2 w l t \"\"\n";
  cout << f << "e\n";
}

template <typename Kernel>
void plot(int up, int down, int& window) {
  Kernel k;
  // Trick to see impulse response
  vector<float> f{0.f, 0.f, 0.f, 0.f, 1.f, .0f, 0.f, 0.f, 0.f};
  const int w = f.size();
  f = upsample(f, up, k);
  gnuplot(f, window++, k, w);
  f = downsample(f, down, k);
  gnuplot(f, window++, k, w);
}

template<typename Kernel>
static void check_interpolation() {
  Kernel k;
  cout << "checking " << k.name() << '\n';
  vector<float> f{.5f, 2.f, 1.f, 0.f, 5.f};
  vector<float> g = upsample(f, f.size(), k);
  cout << "  " << maxerr(f, g) << "  max error\n";
  g = downsample(f, f.size(), k);
  cout << "  " << maxerr(f, g) << "  max error\n";
}


} // end of namespace GS
// int main() {
//   if (1) {
//     const int up = 1000001, down = 101;
//     test_performance<Box>(up, down);
//     test_performance<Hat>(up, down);
//     test_performance<CatmullRom>(up, down);
//     test_performance<MitchellNetravali>(up, down);
//     test_performance<CardinalBspline3>(up, down);
//     test_performance<CardinalOMOMS3>(up, down);
//   } else if (0) {
//     const int up = 1001, down = 51;
//     int window = 1;
//     plot<Box>(up, down, window);
//     plot<Hat>(up, down, window);
//     plot<CatmullRom>(up, down, window);
//     plot<MitchellNetravali>(up, down, window);
//     plot<CardinalBspline3>(up, down, window);
//     plot<CardinalOMOMS3>(up, down, window);
//   } else {
//     check_interpolation<Box>();
//     check_interpolation<Hat>();
//     check_interpolation<CatmullRom>();
//     check_interpolation<CardinalBspline3>();
//     check_interpolation<CardinalOMOMS3>();
//   }
// }