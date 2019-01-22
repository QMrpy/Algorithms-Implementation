#include <bits/stdc++.h>
#define pi 3.14159265
using namespace std;
vector<complex<double> > fft(vector<complex<double> > u,int n)
{
int c,k;
double re,img;
complex<double>cpx(1,0);
complex<double>w(cos(2*pi/n),sin(2*pi/n));
vector<complex<double> >even,odd,even_fft,odd_fft,y;
vector<complex<double> >::iterator j,j1,j2;
if(n==1)
return u;
else
{
c=0;
for(j=u.begin();j!=u.end();j++)
{
if(c%2==0)
even.push_back(complex<double>(real(*j),imag(*j)));
else
odd.push_back((complex<double>(real(*j),imag(*j))));
c++;
}
even_fft=fft(even,n/2);
odd_fft=fft(odd,n/2);
j1=even_fft.begin();j2=odd_fft.begin();
for(k=1;k<=n/2;k++)
{
re=real(*(j1)+cpx*(*(j2)));
img=imag(*(j1)+cpx*(*(j2)));
y.push_back(complex<double>(re,img));
j1++;j2++;
cpx*=w;
}
complex<double>cx(1,0);
j1=even_fft.begin();j2=odd_fft.begin();
for(k=1;k<=n/2;k++)
{
re=real(*(j1)-cx*(*(j2)));
img=imag(*(j1)-cx*(*(j2)));
y.push_back(complex<double>(re,img));
j1++;j2++;
cx*=w;
}
return y;
}
}

vector<complex<double> > inv_fft(vector<complex<double> > u,int n)
{
int c,k;
double re,img;
complex<double>cpx(1,0);
complex<double>w(cos(2*pi/n),-sin(2*pi/n));
vector<complex<double> >even,odd,even_fft,odd_fft,y,z;
vector<complex<double> >::iterator j,j1,j2;
if(n==1)
return u;
else
{
c=0;
for(j=u.begin();j!=u.end();j++)
{
if(c%2==0)
even.push_back(complex<double>(real(*j),imag(*j)));
else
odd.push_back(complex<double>(real(*j),imag(*j)));
c++;
}
even_fft=inv_fft(even,n/2);
odd_fft=inv_fft(odd,n/2);
j1=even_fft.begin();j2=odd_fft.begin();
for(k=1;k<=n/2;k++)
{
re=real(*(j1)+cpx*(*(j2)));
img=imag(*(j1)+cpx*(*(j2)));
y.push_back(complex<double>(re,img));
j1++;j2++;
cpx*=w;
}
complex<double>cx(1,0);
j1=even_fft.begin();j2=odd_fft.begin();
for(k=1;k<=n/2;k++)
{
re=real(*(j1)-cx*(*(j2)));
img=imag(*(j1)-cx*(*(j2)));
y.push_back(complex<double>(re,img));
j1++;j2++;
cx*=w;
}
for(j=y.begin();j!=y.end();j++)
{
re=real(*j);
img=imag(*j);
re=re/n;
img=img/n;
z.push_back(complex<double>(re,img));
}
return z;
}
}

int main()
{
int n,i,a;
vector<complex<double> >v1,v2;
cout<<"Enter polynomial degree as power of 2\n";
cin>>n;
cout<<"Enter polynomial 1\n";
for(i=0;i<=n;i++)
{
cin>>a;
v1.push_back(complex<double>(a,0));
}
vector<complex<double> >p=fft(v1,n);
cout<<"Enter polynomial 2\n";
for(i=0;i<n;i++)
{
cin>>a;
v2.push_back(complex<double>(a,0));
}
vector<complex<double> >q=fft(v2,n);
vector<complex<double> >r;
vector<complex<double> >::iterator j,l;
l=q.begin();
for(j=p.begin();j!=p.end();j++)
{
r.push_back(complex<double>(real((*(j))*(*(l))),imag((*(j))*(*(l)))));
l++;
}
vector<complex<double> >ans=inv_fft(r,2*n);
cout<<"The product of the polynomials is\n";
for(j=ans.begin();j!=ans.end();j++)
cout<<*j<<" ";
cout<<"\n";
return 0;
}



