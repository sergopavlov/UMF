using System;
using System.Collections.Generic;
using System.IO;

namespace UMF3
{
    public enum edges
    {
        Top,
        Bottom,
        Left,
        Right,
        Front,
        Back
    }
    internal class Program
    {
        static void Main(string[] args)
        {
            Func<double, double, double, double> ressin = (x, y, z) => y;
            Func<double, double, double, double> rescos = (x, y, z) => y;
            List<BC> BCs = new();
            for (int i = 0; i < 5; i++)
            {
                BCs.Add(new BC1((edges)i, rescos, 0));
                BCs.Add(new BC1((edges)i, ressin, 1));
            }
            BCs.Add(new BC2(edges.Back, (x, y, z) => 1, 0));
            BCs.Add(new BC2(edges.Back, (x, y, z) => 1, 1));
            Mke mke = new Mke((x, y, z) => 1, (x, y, z) => 1, (x, y, z) => 1, (x, y, z) => 0, (x, y, z) => -2 * y, 1, BCs);
            mke.ReadMesh();
            mke.SolveLU();
            Console.WriteLine(mke.GetSollution(1, 1, 2));
        }
    }
    public class Mke
    {
        int n;
        int m;
        private List<double> Xgrid = new();
        private List<double> Ygrid = new();
        private List<double> Zgrid = new();
        MatrixSparce mat = new();
        MatrixProfile matprof;
        List<double> q = new();
        private Func<double, double, double, double> lambda;
        private Func<double, double, double, double> sigma;
        private Func<double, double, double, double> hi;
        private Func<double, double, double, double> fc;
        private Func<double, double, double, double> fs;
        public List<Elem> elems = new();
        private List<BC> BCs;
        private double omega;
        private int GetGlobalNum(Elem el, int localnum)
        {
            return 2 * el.xlownum + 2 * n * el.ylownum + 2 * m * el.zlownum + localnum % 2
                + (localnum / 2) % 2 * 2 + (localnum / 4) % 2 * n * 2 + localnum / 8 * m * 2;
        }
        public Mke(Func<double, double, double, double> lambda, Func<double, double, double, double> sigma, Func<double, double, double, double> hi, Func<double, double, double, double> fc, Func<double, double, double, double> fs, double omega, List<BC> BCs)
        {
            this.lambda = lambda;
            this.sigma = sigma;
            this.hi = hi;
            this.fc = fc;
            this.fs = fs;
            this.omega = omega;
            this.BCs = BCs;
        }
        public void ReadMesh()
        {
            foreach (var item in File.ReadAllLines("TXT/Xgrid.txt"))
            {
                Xgrid.Add(double.Parse(item));
            }
            foreach (var item in File.ReadAllLines("TXT/Ygrid.txt"))
            {
                Ygrid.Add(double.Parse(item));
            }
            foreach (var item in File.ReadAllLines("TXT/Zgrid.txt"))
            {
                Zgrid.Add(double.Parse(item));
            }
            n = Xgrid.Count;
            m = Xgrid.Count * Ygrid.Count;
            for (int i = 0; i < n - 1; i++)
            {
                for (int j = 0; j < Ygrid.Count - 1; j++)
                {
                    for (int p = 0; p < Zgrid.Count - 1; p++)
                    {
                        elems.Add(new Elem(i, j, p));
                    }
                }
            }
        }
        public void GenerateProfile()
        {
            mat.n = 2 * Xgrid.Count * Ygrid.Count * Zgrid.Count;
            mat.di = new List<double>(mat.n);
            mat.ia = new List<int>(mat.n + 1);
            mat.ja = new List<int>();
            mat.al = new List<double>();
            mat.au = new List<double>();
            mat.b = new List<double>(mat.n);
            List<SortedSet<int>> list = new List<SortedSet<int>>(mat.n);
            for (int i = 0; i < mat.n; i++)
            {
                list.Add(new SortedSet<int>());
            }
            foreach (var item in elems)
            {
                for (int i = 0; i < 16; i++)
                {
                    for (int j = 0; j < 16; j++)
                    {
                        list[GetGlobalNum(item, i)].Add(GetGlobalNum(item, j));
                    }
                }
            }
            mat.ia.Add(0);
            mat.ia.Add(0);
            mat.di.Add(0);
            mat.b.Add(0);
            for (int i = 1; i < mat.n; i++)
            {
                mat.di.Add(0);
                mat.b.Add(0);
                int count = 0;
                foreach (var item in list[i])
                {
                    if (item < i)
                    {
                        mat.ja.Add(item);
                        mat.al.Add(0);
                        mat.au.Add(0);
                        count++;
                    }
                }
                mat.ia.Add(mat.ia[i] + count);
            }
        }
        public void Addlocal(Elem el)
        {
            double hx = Xgrid[el.xlownum + 1] - Xgrid[el.xlownum];
            double hy = Ygrid[el.ylownum + 1] - Ygrid[el.ylownum];
            double hz = Zgrid[el.zlownum + 1] - Zgrid[el.zlownum];
            double Xmean = (Xgrid[el.xlownum + 1] + Xgrid[el.xlownum]) / 2;
            double Ymean = (Ygrid[el.ylownum + 1] + Ygrid[el.ylownum]) / 2;
            double Zmean = (Zgrid[el.zlownum + 1] + Zgrid[el.zlownum]) / 2;
            double[] fcos = new double[8];
            double[] fsin = new double[8];
            for (int i = 0; i < 8; i++)
            {
                fcos[i] = fc(Xgrid[el.xlownum + i % 2], Ygrid[el.ylownum + (i / 2) % 2], Zgrid[el.zlownum + i / 4]);
                fsin[i] = fs(Xgrid[el.xlownum + i % 2], Ygrid[el.ylownum + (i / 2) % 2], Zgrid[el.zlownum + i / 4]);
            }
            for (int i = 0; i < 16; i += 2)//заполняем косинусовые строчки
            {
                double curb = 0;
                for (int j = 0; j < 8; j++)
                {
                    curb += fcos[j] * Matrices.Mij(i / 2, j, hx, hy, hz);
                }
                mat.b[GetGlobalNum(el, i)] += curb;
                mat.di[GetGlobalNum(el, i)] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, i / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, i / 2, hx, hy, hz);
                for (int j = 0; j < i; j += 2)//заполняем cos cos
                {
                    int max = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, i) : GetGlobalNum(el, j);
                    int min = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, j) : GetGlobalNum(el, i);
                    int k = mat.ia[max];
                    while (mat.ja[k] < min)
                    {
                        k++;
                    }
                    mat.al[k] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, j / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                    mat.au[k] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, j / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                }
                for (int j = 1; j < i; j += 2)//заполняем cos sin
                {
                    int max = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, i) : GetGlobalNum(el, j);
                    int min = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, j) : GetGlobalNum(el, i);
                    int k = mat.ia[max];
                    while (mat.ja[k] < min)
                    {
                        k++;
                    }
                    mat.al[k] += omega * sigma(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                    mat.au[k] += -omega * sigma(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                }
            }
            for (int i = 1; i < 16; i += 2)//заполняем синусовые строчки
            {
                double curb = 0;
                for (int j = 0; j < 8; j++)
                {
                    curb += fsin[j] * Matrices.Mij(i / 2, j, hx, hy, hz);
                }
                mat.b[GetGlobalNum(el, i)] += curb;
                mat.di[GetGlobalNum(el, i)] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, i / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, i / 2, hx, hy, hz);
                for (int j = 0; j < i; j += 2)// заполняем sin cos
                {
                    int max = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, i) : GetGlobalNum(el, j);
                    int min = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, j) : GetGlobalNum(el, i);
                    int k = mat.ia[max];
                    while (mat.ja[k] < min)
                    {
                        k++;
                    }
                    mat.al[k] += -omega * sigma(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                    mat.au[k] += omega * sigma(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                }
                for (int j = 1; j < i; j += 2)// заполняем sin sin
                {
                    int max = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, i) : GetGlobalNum(el, j);
                    int min = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, j) : GetGlobalNum(el, i);
                    int k = mat.ia[max];
                    while (mat.ja[k] < min)
                    {
                        k++;
                    }
                    mat.al[k] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, j / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                    mat.au[k] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, j / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                }
            }
        }
        public void AddBoundary3()
        {
            foreach (var item in BCs.FindAll((T) => T.bctype == 3))
            {
                var curitem = item as BC3;
                switch (curitem.edge)
                {
                    case edges.Top:
                        foreach (var el in elems.FindAll((T) => T.zlownum + 1 == Zgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 8), GetGlobalNum(el, 10), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 9), GetGlobalNum(el, 11), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[Zgrid.Count - 1];
                            double hx = x[1] - x[0], hy = y[1] - y[0];

                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x[j % 2], y[j / 2], z) * Matrices.Mij2(i, j, hx, hy);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hx, hy);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hx, hy);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hx, hy);
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Bottom:
                        foreach (var el in elems.FindAll((T) => T.zlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 4), GetGlobalNum(el, 6) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 5), GetGlobalNum(el, 7) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[0];
                            double hx = x[1] - x[0], hy = y[1] - y[0];

                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x[j % 2], y[j / 2], z) * Matrices.Mij2(i, j, hx, hy);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hx, hy);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hx, hy);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hx, hy);
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Left:
                        foreach (var el in elems.FindAll((T) => T.xlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 4), GetGlobalNum(el, 8), GetGlobalNum(el, 12) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 5), GetGlobalNum(el, 9), GetGlobalNum(el, 13) };
                            }
                            double x = Xgrid[0];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hy = y[1] - y[0];

                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x, y[j % 2], z[j / 2]) * Matrices.Mij2(i, j, hz, hy);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hz, hy);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hz, hy);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hz, hy);
                                    }
                                }
                            }


                        }
                        break;
                    case edges.Right:
                        foreach (var el in elems.FindAll((T) => T.xlownum + 1 == Xgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 2), GetGlobalNum(el, 6), GetGlobalNum(el, 10), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 3), GetGlobalNum(el, 7), GetGlobalNum(el, 11), GetGlobalNum(el, 15) };
                            }
                            double x = Xgrid[Xgrid.Count - 1];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hy = y[1] - y[0];
                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x, y[j % 2], z[j / 2]) * Matrices.Mij2(i, j, hz, hy);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hz, hy);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hz, hy);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hz, hy);
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Front:
                        foreach (var el in elems.FindAll((T) => T.ylownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 8), GetGlobalNum(el, 10) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 9), GetGlobalNum(el, 11) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[0];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hx = x[1] - x[0];

                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x[j % 2], y, z[j / 2]) * Matrices.Mij2(i, j, hz, hz);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hz, hz);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hz, hz);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hz, hz);
                                    }
                                }
                            }

                        }
                        break;
                    case edges.Back:
                        foreach (var el in elems.FindAll((T) => T.ylownum + 1 == Ygrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 4), GetGlobalNum(el, 6), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 5), GetGlobalNum(el, 7), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[Ygrid.Count - 1];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hx = x[1] - x[0];
                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x[j % 2], y, z[j / 2]) * Matrices.Mij2(i, j, hz, hz);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hz, hz);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hz, hz);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hz, hz);
                                    }
                                }
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        public void AddBoundary2()
        {
            foreach (var item in BCs.FindAll((T) => T.bctype == 2))
            {
                var curitem = item as BC2;
                switch (curitem.edge)
                {
                    case edges.Top:
                        foreach (var el in elems.FindAll((T) => T.zlownum + 1 == Zgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 8), GetGlobalNum(el, 10), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 9), GetGlobalNum(el, 11), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[Zgrid.Count - 1];
                            double hx = x[1] - x[0], hy = y[1] - y[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x[j % 2], y[j / 2], z) * Matrices.Mij2(i, j, hx, hy);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    case edges.Bottom:
                        foreach (var el in elems.FindAll((T) => T.zlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 4), GetGlobalNum(el, 6) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 5), GetGlobalNum(el, 7) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[0];
                            double hx = x[1] - x[0], hy = y[1] - y[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x[j % 2], y[j / 2], z) * Matrices.Mij2(i, j, hx, hy);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    case edges.Left:
                        foreach (var el in elems.FindAll((T) => T.xlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 4), GetGlobalNum(el, 8), GetGlobalNum(el, 12) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 5), GetGlobalNum(el, 9), GetGlobalNum(el, 13) };
                            }
                            double x = Xgrid[0];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hy = y[1] - y[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x, y[j % 2], z[j / 2]) * Matrices.Mij2(i, j, hz, hy);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    case edges.Right:
                        foreach (var el in elems.FindAll((T) => T.xlownum + 1 == Xgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 2), GetGlobalNum(el, 6), GetGlobalNum(el, 10), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 3), GetGlobalNum(el, 7), GetGlobalNum(el, 11), GetGlobalNum(el, 15) };
                            }
                            double x = Xgrid[Xgrid.Count - 1];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hy = y[1] - y[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x, y[j % 2], z[j / 2]) * Matrices.Mij2(i, j, hz, hy);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    case edges.Front:
                        foreach (var el in elems.FindAll((T) => T.ylownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 8), GetGlobalNum(el, 10) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 9), GetGlobalNum(el, 11) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[0];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hx = x[1] - x[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x[j % 2], y, z[j / 2]) * Matrices.Mij2(i, j, hz, hx);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    case edges.Back:
                        foreach (var el in elems.FindAll((T) => T.ylownum + 1 == Ygrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 4), GetGlobalNum(el, 6), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 5), GetGlobalNum(el, 7), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[Ygrid.Count - 1];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hx = x[1] - x[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x[j % 2], y, z[j / 2]) * Matrices.Mij2(i, j, hz, hx);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        public void AddBoundary1()
        {
            foreach (var item in BCs.FindAll((T) => T.bctype == 1))
            {
                var curitem = item as BC1;
                switch (curitem.edge)
                {
                    case edges.Top:
                        foreach (var el in elems.FindAll((T) => T.zlownum + 1 == Zgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 8), GetGlobalNum(el, 10), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 9), GetGlobalNum(el, 11), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[Zgrid.Count - 1];
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x[p % 2], y[p / 2], z);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x[p % 2], y[p / 2], z) * mat.au[k];
                                    mat.al[k] = 0;
                                    mat.au[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int k = mat.ia[i];
                                    while (mat.ja[k] < index[p] && k < mat.ia[i + 1] - 1)
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[p])
                                    {
                                        mat.b[i] -= mat.al[k] * curitem.value(x[p % 2], y[p / 2], z);
                                        mat.au[k] = 0;
                                        mat.al[k] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Bottom:
                        foreach (var el in elems.FindAll((T) => T.zlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 4), GetGlobalNum(el, 6) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 5), GetGlobalNum(el, 7) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[0];
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x[p % 2], y[p / 2], z);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x[p % 2], y[p / 2], z) * mat.au[k];
                                    mat.au[k] = 0;
                                    mat.al[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int k = mat.ia[i];
                                    while (mat.ja[k] < index[p] && k < mat.ia[i + 1] - 1)
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[p])
                                    {
                                        mat.b[i] -= mat.al[k] * curitem.value(x[p % 2], y[p / 2], z);
                                        mat.al[k] = 0;
                                        mat.au[k] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Left:
                        foreach (var el in elems.FindAll((T) => T.xlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 4), GetGlobalNum(el, 8), GetGlobalNum(el, 12) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 5), GetGlobalNum(el, 9), GetGlobalNum(el, 13) };
                            }
                            double x = Xgrid[0];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x, y[p % 2], z[p / 2]);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x, y[p % 2], z[p / 2]) * mat.au[k];
                                    mat.al[k] = 0;
                                    mat.au[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int k = mat.ia[i];
                                    while (mat.ja[k] < index[p] && k < mat.ia[i + 1] - 1)
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[p])
                                    {
                                        mat.b[i] -= mat.al[k] * curitem.value(x, y[p % 2], z[p / 2]);
                                        mat.au[k] = 0;
                                        mat.al[k] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Right:
                        foreach (var el in elems.FindAll((T) => T.xlownum + 1 == Xgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 2), GetGlobalNum(el, 6), GetGlobalNum(el, 10), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 3), GetGlobalNum(el, 7), GetGlobalNum(el, 11), GetGlobalNum(el, 15) };
                            }
                            double x = Xgrid[Xgrid.Count - 1];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x, y[p % 2], z[p / 2]);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x, y[p % 2], z[p / 2]) * mat.au[k];
                                    mat.al[k] = 0;
                                    mat.au[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int k = mat.ia[i];
                                    while (mat.ja[k] < index[p] && k < mat.ia[i + 1] - 1)
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[p])
                                    {
                                        mat.b[i] -= mat.al[k] * curitem.value(x, y[p % 2], z[p / 2]);
                                        mat.au[k] = 0;
                                        mat.al[k] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Front:
                        foreach (var el in elems.FindAll((T) => T.ylownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 8), GetGlobalNum(el, 10) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 9), GetGlobalNum(el, 11) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[0];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x[p % 2], y, z[p / 2]);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x[p % 2], y, z[p / 2]) * mat.au[k];
                                    mat.al[k] = 0;
                                    mat.au[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int k = mat.ia[i];
                                    while (mat.ja[k] < index[p] && k < mat.ia[i + 1] - 1)
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[p])
                                    {
                                        mat.b[i] -= mat.al[k] * curitem.value(x[p % 2], y, z[p / 2]);
                                        mat.au[k] = 0;
                                        mat.al[k] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Back:
                        foreach (var el in elems.FindAll((T) => T.ylownum + 1 == Ygrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 4), GetGlobalNum(el, 6), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 5), GetGlobalNum(el, 7), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[Ygrid.Count - 1];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x[p % 2], y, z[p / 2]);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x[p % 2], y, z[p / 2]) * mat.au[k];
                                    mat.al[k] = 0;
                                    mat.au[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int k = mat.ia[i];
                                    while (mat.ja[k] < index[p] && k < mat.ia[i + 1] - 1)
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[p])
                                    {
                                        mat.b[i] -= mat.al[k] * curitem.value(x[p % 2], y, z[p / 2]);
                                        mat.au[k] = 0;
                                        mat.al[k] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        public void SolveLU()
        {
            GenerateProfile();
            foreach (var item in elems)
            {
                Addlocal(item);
            }
            AddBoundary2();
            AddBoundary3();
            AddBoundary1();
            matprof = mat.ToProfile();
            q = matprof.Solve();
        }
        public void SolveLOS(double eps, int maxiter)
        {
            GenerateProfile();
            foreach (var item in elems)
            {
                Addlocal(item);
            }
            AddBoundary2();
            AddBoundary3();
            AddBoundary1();
            List<double> x0 = new();
            for (int i = 0; i < mat.n; i++)
            {
                x0.Add(0);
            }
            q = mat.LOS(x0, maxiter, eps);
        }
        public (double, double) GetSollution(double x, double y, double z)
        {
            bool flag = true;
            int n = elems.Count;
            int i = 0;
            for (i = 0; i < n && flag; i++)
            {
                if (x >= Xgrid[elems[i].xlownum] && x <= Xgrid[elems[i].xlownum + 1] && y >= Ygrid[elems[i].ylownum] && y <= Ygrid[elems[i].ylownum + 1] && z >= Zgrid[elems[i].zlownum] && z <= Zgrid[elems[i].zlownum + 1])
                {
                    flag = false;
                    i--;
                }
            }
            if (flag)
                throw new Exception();
            double xloc = (x - Xgrid[elems[i].xlownum]) / (Xgrid[elems[i].xlownum + 1] - Xgrid[elems[i].xlownum]);
            double yloc = (y - Ygrid[elems[i].ylownum]) / (Ygrid[elems[i].ylownum + 1] - Ygrid[elems[i].ylownum]);
            double zloc = (z - Zgrid[elems[i].zlownum]) / (Zgrid[elems[i].zlownum + 1] - Zgrid[elems[i].zlownum]);
            double[] psix = { 1 - xloc, xloc };
            double[] psiy = { 1 - yloc, yloc };
            double[] psiz = { 1 - zloc, zloc };
            (double, double) res = (0, 0);
            for (int p = 0; p < 8; p++)
            {
                res.Item1 += q[GetGlobalNum(elems[i], 2 * p)] * psix[p % 2] * psiy[(p / 2) % 2] * psiz[p / 4];
                res.Item2 += q[GetGlobalNum(elems[i], 2 * p + 1)] * psix[p % 2] * psiy[(p / 2) % 2] * psiz[p / 4];
            }
            return res;
        }
    }
    public class Elem
    {
        public Elem(int xlownum, int ylownum, int zlownum)
        {
            this.xlownum = xlownum;
            this.ylownum = ylownum;
            this.zlownum = zlownum;
        }

        public int xlownum { get; init; }
        public int ylownum { get; init; }
        public int zlownum { get; init; }
    }
    public abstract class BC
    {
        public edges edge { get; init; }
        public int cossintype { get; init; }
        public int bctype { get; init; }

        public Func<double, double, double, double> value;

        public BC(int bctype, edges edge, Func<double, double, double, double> value, int cossintype)
        {
            this.bctype = bctype;
            this.edge = edge;
            this.value = value;
            this.cossintype = cossintype;
        }
    }
    public class BC1 : BC
    {
        public BC1(edges edge, Func<double, double, double, double> value, int cossintype) : base(1, edge, value, cossintype) { }
    }
    public class BC2 : BC
    {
        public BC2(edges edge, Func<double, double, double, double> value, int cossintype) : base(2, edge, value, cossintype) { }
    }
    public class BC3 : BC
    {
        public double betta { get; init; }
        public BC3(edges edge, Func<double, double, double, double> value, double betta, int cossintype) : base(3, edge, value, cossintype) { this.betta = betta; }
    }
    public static class Matrices
    {
        static double[][] Gmatr = new double[][]
        {
            new double[]
            {
                1,-1
            },

            new double[]
            {
                -1,1
            }
        };
        static double[][] Mmatr = new double[][]
       {
            new double[]
            {
               2.0/6,1.0/6
            },
            new double[]
            {
                1.0/6,2.0/6
            }
       };
        public static double Gij(int i, int j, double hx, double hy, double hz)
        {
            return hy * hz / hx * Gmatr[i % 2][j % 2] * Mmatr[(i / 2) % 2][(j / 2) % 2] * Mmatr[i / 4][j / 4] + hx * hz / hy * Mmatr[i % 2][j % 2] * Gmatr[(i / 2) % 2][(j / 2) % 2] * Mmatr[i / 4][j / 4] + hx * hy / hz * Mmatr[i % 2][j % 2] * Mmatr[(i / 2) % 2][(j / 2) % 2] * Gmatr[i / 4][j / 4];
        }
        public static double Mij(int i, int j, double hx, double hy, double hz) =>
            hx * hy * hz * Mmatr[i % 2][j % 2] * Mmatr[(i / 2) % 2][(j / 2) % 2] * Mmatr[i / 4][j / 4];
        public static double Mij2(int i, int j, double hx, double hy) =>
            hx * hy * Mmatr[i % 2][j % 2] * Mmatr[i / 2][j / 2];
    }
    public class Matrix
    {

    }
    public class MatrixSparce
    {
        public List<double> al;
        public List<double> au;
        public List<double> di;
        public List<int> ia;
        public List<int> ja;
        public List<double> b;
        public int n;
        public MatrixSparce() { }
        public List<double> MSG(List<double> x0, int maxiter, double eps)
        {
            double bnorm = Math.Sqrt(DotProduct(b, b));
            List<double> p = new List<double>(n);
            List<double> q = new List<double>(n);
            for (int i = 0; i < n; i++)
            {
                p.Add(0);
            }
            var r = MatrixMult(x0);
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - r[i];
                p[i] = r[i];
            }
            int k = 0;
            double alpha, betta, rnorm = Math.Sqrt(DotProduct(r, r));
            while (k < maxiter && rnorm / bnorm > eps)
            {
                q = MatrixMult(p);
                alpha = DotProduct(r, r) / DotProduct(q, p);
                betta = 1 / DotProduct(r, r);
                for (int i = 0; i < n; i++)
                {
                    x0[i] += alpha * p[i];
                    r[i] -= alpha * q[i];
                }
                rnorm = Math.Sqrt(DotProduct(r, r));
                betta *= DotProduct(r, r);
                for (int i = 0; i < n; i++)
                {
                    p[i] = r[i] + betta * p[i];
                }
                k++;
            }
            return x0;
        }
        public List<double> LOS(List<double> x0, int maxiter, double eps)
        {
            double bnorm = Math.Sqrt(DotProduct(b, b));

            List<double> Ar = new();
            List<double> r = MatrixMult(x0);
            List<double> z = new();
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - r[i];
                z.Add(r[i]);
            }
            List<double> p = MatrixMult(z);
            int k = 0;
            double alpha, betta, rnorm = Math.Sqrt(DotProduct(r, r));
            while (k < maxiter && rnorm / bnorm > eps)
            {
                alpha = DotProduct(p, r) / DotProduct(p, p);
                for (int i = 0; i < n; i++)
                {
                    x0[i] += alpha * z[i];
                    r[i] -= alpha * p[i];
                }
                Ar = MatrixMult(r);
                betta = -DotProduct(p, Ar) / DotProduct(p, p);
                rnorm = Math.Sqrt(DotProduct(r, r));
                for (int i = 0; i < n; i++)
                {
                    z[i] = r[i] + betta * z[i];
                    p[i] = Ar[i] + betta * p[i];
                }
                k++;
            }
            Console.WriteLine($"{rnorm / bnorm} {k}");
            return x0;
        }
        public void PrintDense()
        {
            for (int i = 0; i < n; i++)
            {
                int k = ia[i];
                for (int j = 0; j < i; j++)
                {
                    while (ja[k] < j && k < ia[i + 1] - 1)
                        k++;
                    if (ja[k] == j)
                        Console.Write($"{al[k]:f2} ");
                    else
                        Console.Write($"{0:f2} ");
                }
                Console.Write($"{di[i]:f2} ");
                for (int j = i + 1; j < n; j++)
                {
                    k = ia[j];
                    while (ja[k] < i && k < ia[i + 1] - 1)
                        k++;
                    if (ja[k] == i)
                        Console.Write($"{au[k]:f2} ");
                    else
                        Console.Write($"{0:f2} ");
                }

                Console.WriteLine();
            }
        }
        public MatrixProfile ToProfile()
        {
            MatrixProfile res = new();
            res.n = n;
            res.di = di;
            res.b = b;
            res.ia = new();
            res.al = new();
            res.au = new();
            res.ia.Add(0);
            res.ia.Add(0);
            for (int i = 1; i < n; i++)
            {
                if (ia[i + 1] - ia[i] > 0)
                {
                    int count = 0;
                    int k = ia[i];
                    for (int p = ja[ia[i]]; p < i; p++)
                    {
                        count++;
                        while (ja[k] < p)
                            k++;
                        if (ja[k] == p)
                        {
                            res.al.Add(al[k]);
                            res.au.Add(au[k]);
                        }
                        else
                        {
                            res.al.Add(0);
                            res.au.Add(0);
                        }
                    }
                    res.ia.Add(res.ia[res.ia.Count - 1] + count);
                }
                else
                {
                    res.ia.Add(res.ia[i - 1]);
                }
            }
            return res;
        }
        double DotProduct(List<double> vec1, List<double> vec2)
        {
            if (vec1.Count != vec2.Count)
                throw new Exception();
            double res = 0;
            int m = vec1.Count;
            for (int i = 0; i < m; i++)
            {
                res += vec1[i] * vec2[i];
            }
            return res;
        }
        List<double> MatrixMult(List<double> x)
        {
            if (x.Count != n)
                throw new Exception();
            List<double> res = new List<double>(x.Count);
            for (int i = 0; i < n; i++)
            {
                res.Add(0);
            }
            for (int i = 0; i < n; i++)
            {
                res[i] = x[i] * di[i];
                for (int k = ia[i]; k < ia[i + 1]; k++)
                {
                    int j = ja[k];
                    res[i] += al[k] * x[j];
                    res[j] += au[k] * x[i];
                }
            }
            return res;
        }
    }
    public class MatrixProfile
    {
        public List<double> al;
        public List<double> au;
        public List<double> di;
        public List<int> ia;
        public List<double> b;
        public int n;
        public MatrixProfile() { }
        void LU()
        {
            for (int i = 0; i < n; i++)
            {
                double sumrow = 0;
                int j0 = i - (ia[i + 1] - ia[i]);
                int dind = j0;
                for (int ii = ia[i]; ii < ia[i + 1]; ii++, dind++)
                {
                    int j = j0 + ii - ia[i];
                    int j0j = j - (ia[j + 1] - ia[j]);

                    int kbeg = j0 > j0j ? j0 : j0j;
                    int kend = i > j ? j : i;
                    double summu = 0;
                    double summl = 0;
                    int indexi = ia[i] + kbeg - j0;
                    int indexj = ia[j] + kbeg - j0j;
                    for (int k = 0; k < kend - kbeg; k++, indexi++, indexj++)
                    {
                        summu += al[indexj] * au[indexi];
                        summl += al[indexi] * au[indexj];
                    }
                    al[ii] -= summl;
                    au[ii] = (au[ii] - summu) / di[j];
                    sumrow += au[ii] * al[ii];
                }
                di[i] -= sumrow;
            }
        }
        void GaussU()
        {
            for (int i = n - 1; i >= 0; i--)
            {
                double cursum = 0;
                for (int j = n - 1; j > i; j--)
                {
                    if (ia[j + 1] - ia[j] >= j - i)
                    {
                        cursum += au[ia[j + 1] - j + i] * b[j];
                    }
                }
                b[i] -= cursum;
            }
        }
        void GaussL()
        {
            for (int i = 0; i < n; i++)
            {
                double cursum = 0;
                int j = i - (ia[i + 1] - ia[i]);
                for (int k = ia[i]; k < ia[i + 1]; k++, j++)
                {
                    cursum += al[k] * b[j];
                }
                b[i] = (b[i] - cursum) / di[i];
            }
        }
        public List<double> Solve()
        {
            LU();
            GaussL();
            GaussU();
            return b;
        }
    }

}
