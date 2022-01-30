#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#define N (2<<18)
#define eps 1e-7
using namespace std;
int tint(double x)
{
    int k=(int)(x);
    //cout<<x<<' '<<k<<endl;
    if (abs(x-k)>eps)
    if (k>0||k==0&&x>0) ++k;
    else {
        --k;
    }
    return k;
}
const double pi=acos(-1);

class Complex
{
public:
    double r,i;
    Complex(double _r,double _i):r(_r),i(_i){}
    Complex(){}
};
ostream& operator<<(ostream &output,Complex a)
{
    output<<tint(a.r)<<'+'<<tint(a.i)<<'i';
}
Complex b[N];
inline Complex operator+(const Complex &a,const Complex &b) { return Complex(a.r+b.r,a.i+b.i);}
inline Complex operator*(const Complex &a,const Complex &b) { return Complex(a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r);}
inline Complex operator-(const Complex &a,const Complex &b) { return Complex(a.r-b.r,a.i-b.i);}
inline Complex add(const Complex &a,const Complex &b) { return Complex(a.r+b.r,a.i+b.i);}
inline Complex mul(const Complex &a,const Complex &b) { return Complex(a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r);}
inline Complex sub(const Complex &a,const Complex &b) { return Complex(a.r-b.r,a.i-b.i);}
Complex div(const Complex &a,const Complex &b) { return Complex((a.r*b.r+a.i*b.i)/(b.r*b.r+b.i*b.i),(a.i*b.r-a.r*b.i)/(b.r*b.r+b.i*b.i));}
Complex ome(int i,int n) { return Complex(cos(i*2*pi/n),sin(i*2*pi/n));}
int Rev[N];
int rev(int x,int k)
{
    int r=0,b;
    for (int i=0;i<k;++i)
    {
        b=(x>>i)&1;
        r+=(b<<(k-i-1));
    }
    return r;
}
Complex w[2][N];
int z=1;
void dft(Complex *a,int k,int type)
{
    for (int i=0;i<(1<<k);++i) b[i]=a[Rev[i]];
    //for (int j=0;j<(1<<k);++j) cout<<b[j]<<' ';
    //cout<<endl;
    for (int i=1;i<=k;++i)
    {
        //for (int j=0;j<(1<<k);++j) { a[j]=add(b[j&(~(1<<(i-1)))],mul(ome(type*j,(1<<i)),b[j|(1<<(i-1))]));}
        //for (int j=0;j<(1<<k);++j) b[j]=a[j];
        for (int j=0;j<(1<<k);++j)
        {
            a[j]=b[j&(~(1<<(i-1)))]+w[type][(j%(1<<i))*(z/(1<<i))]*b[j|(1<<(i-1))];
            //cout<<w[type][(j%(1<<i))*(z/(1<<i))]-ome((type?-1:1)*j,1<<i)<<endl;
        }//ome(type*j,(1<<i))*b[j|(1<<(i-1))];}
        for (int j=0;j<(1<<k);++j) b[j]=a[j];
        //for (int j=0;j<(1<<k);++j) cout<<b[j]<<' ';
        //cout<<endl;
    }
}
/*void dft(Complex* a,int k,int f)
{
	for(int i=0;i<z;i++) if(Rev[i]>i) swap(a[i],a[Rev[i]]);

	Complex x,y;
	for(int i=1;i<z;i<<=1) for(int j=0,l=z/(i<<1);j<z;j+=(i<<1)) for(int k=0,t=0;k<i;k++,t+=l)
		x=a[j+k],y=w[f][t]*a[j+k+i],a[j+k]=x+y,a[j+k+i]=x-y;

	//if(f) for(int i=0;i<n;i++) a[i].r/=n;
}*/
void idft(Complex *a,int k)
{
    for (int i=0;i<(1<<k);++i) b[i]=a[i];
    for (int i=k-1;i>=0;--i)
    {
        for (int j=0;j<(1<<k);++j) if (j%(1<<(i+1))<(1<<i))
        {
            a[j]=add(b[j],b[j+(1<<i)]);
            a[j+(1<<i)]=div(sub(b[j],b[j+(1<<i)]),ome((j%(1<<i)),1<<(i+1)));
        }
        for (int j=0;j<(1<<k);++j) b[j]=a[j];
        //cout<<endl;
    }
    for (int i=0;i<(1<<k);++i) a[i]=b[Rev[i]];
}
Complex a1[N],a2[N];
int n,m,k=0;
const int base = 2;

void polynomialToInt(vector<int> &res) {
    int i = 0;
    int x = 0;
    int len = res.size();
    while (i < len || x) {
        if (i < len) {
            x += res[i];
        }
        if (i < len) {
            res[i] = x % base;
        } else {
            res.push_back(x % base);
        }
        x /= base; i += 1;
    }
}
vector<int> FFTMulti(vector<int> &a, vector<int> &b)
{
    n = int(a.size()) - 1; m = int(b.size()) - 1;
    // scanf("%d%d",&n,&m);
    for (;z<n+m+1;z*=2,++k);
    double r;
    //for (int i=0;i<z;++i) { w[0][i]=w[1][i]=ome(i,z); w[1][i].i=-w[0][i].i;}
    for(int i=0;i<z;i++) w[0][i]=w[1][i]=Complex(cos(2*pi*i/z),sin(2*pi*i/z)),w[1][i].i=-w[0][i].i;
    for (int i=0;i<z;++i) { a1[i]=a2[i]=Complex(0,0); Rev[i]=rev(i,k);}
    for (int i=0;i<=n;++i)
    {
        // scanf("%lf",&r);
        a1[i].r = a[i];
    }
    for (int i=0;i<=m;++i)
    {
        // scanf("%lf",&r);
        a2[i].r = b[i];
    }
    dft(a1,k,0);dft(a2,k,0);
    //cout<<endl;
    //for (int i=0;i<z;++i) cout<<a1[i]<<' ';
    //cout<<endl;
    for (int i=0;i<z;++i) a1[i]=mul(a1[i],a2[i]);
    dft(a1,k,1);
    vector<int> res(n + m + 1, 0);
    for (int i = 0; i <= n + m; ++i) {
        res[i] = tint(a1[i].r / z);
    }
    // for (int i=0;i<m+n;++i) printf("%d ",tint(a1[i].r/z));
    // printf("%d\n",tint(a1[m+n].r/z));

    polynomialToInt(res);
    return res;
}

vector<int> naiveMulti(vector<int> &a, vector<int> &b) {
    n = int(a.size()) - 1; m = int(b.size()) - 1;
    vector<int> res(n + m + 1, 0);
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            res[i + j] += a[i] * b[j];
        }
    }
    polynomialToInt(res);
    return res;
}
int main()
{
    //freopen("input.txt","r",stdin);
    //freopen("output.txt","w",stdout);
    //cout<<rev(5,3)<<endl;
    //cout<<tint(-1.6)<<endl;

    vector<int> a, b;
    int n, m; 
    cin >> n;
    a = vector<int>(n + 1);
    for (int i = 0; i <= n; ++i) {
        cin >> a[i];
    }
    cin >> m;
    b = vector<int>(m + 1);
    for (int i = 0; i <= m; ++i) {
        cin >> b[i];
    }

    vector<int> res = FFTMulti(a, b);
    for (int i = 0; i < res.size(); ++i) {
        cout << res[i] << ' ';
    }

    cout << "(FFT multi)" << endl;

    vector<int> naiveRes = naiveMulti(a, b);
    for (int i = 0; i < naiveRes.size(); ++i) {
        cout << naiveRes[i] << ' ';
    }

    cout << "(naive multi)" << endl;

    bool flag = true;
    if (res.size() != naiveRes.size()) {
        flag = false;
    }
    for (int i = 0; i < res.size(); ++i) {
        if (res[i] != naiveRes[i]) {
            flag = false; break;
        }
    }

    cout << (flag ? "correct!" : "wrong!") << endl;


    // vector<int> res = 
    return 0;
}