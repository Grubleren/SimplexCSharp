using System;
using System.Windows.Forms;
using System.IO;
using System.Drawing;
using System.Globalization;

namespace SimplexCSharp
{
    public partial class Form1 : Form
    {
        double epsilon = 0.0000000000001;
        new struct Bounds
        {
            public double x;
            public double y;
        }

        int count;

        int noRows;
        int noColumns;
        int rightSideIndex;
        int noArtificials;

        double[][] tableau;
        double[] zFunction;

        double[] wFunction;

        Random random = new Random(2730);

        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            string[] line;
            string inputDataPath = @"C:..\..\TestData5.txt";
            StreamReader reader = new StreamReader(inputDataPath);

            string s = reader.ReadLine();
            string[] size = s.Split(',');
            int noConstraints = int.Parse(size[0]);
            int noVariables = int.Parse(size[1]);
            noColumns = noVariables + 1;

            Bounds[] boundaries = new Bounds[noVariables];
            s = reader.ReadLine();
            string[] limits = s.Split(',');

            int countUpperLimit = 0;
            for (int i = 0; i < noVariables; i++)
            {
                boundaries[i].x = !limits[2 * i].Contains("inf") ? double.Parse(limits[2 * i], CultureInfo.InvariantCulture) : double.MinValue;
                boundaries[i].y = !limits[2 * i + 1].Contains("inf") ? double.Parse(limits[2 * i + 1], CultureInfo.InvariantCulture) : double.MaxValue;
                if (boundaries[i].y != double.MaxValue)
                    countUpperLimit++;
            }

            noRows = noConstraints + countUpperLimit;

            noArtificials = noRows;
            int noSlags = noRows;
            rightSideIndex = noArtificials + noSlags + noColumns - 1;

            zFunction = new double[noArtificials + noSlags + noVariables + 1];
            wFunction = new double[noArtificials + noSlags + noVariables + 1];
            tableau = new double[noRows][];
            for (int i = 0; i < noRows; i++)
                tableau[i] = new double[noArtificials + noSlags + noVariables + 1];

            double c;
            for ( int i = 0; i < noConstraints; i++)
            {
                s = reader.ReadLine();
                line = s.Split(',');
                c = 0;
                for (int j = 0; j < noVariables; j++)
                {
                    double d = double.Parse(line[j], CultureInfo.InvariantCulture);
                    tableau[i][noArtificials + noSlags + j] = d;
                    c += boundaries[j].x != double.MinValue ? boundaries[j].x * d : 0;
                }
                tableau[i][noArtificials + noSlags + noVariables] = double.Parse(line[noVariables + 1], CultureInfo.InvariantCulture) - c;

                if (line[noVariables].Contains("<"))
                    tableau[i][noArtificials + i] = 1;
                if (line[noVariables].Contains(">"))
                    tableau[i][noArtificials + i] = -1;
                if (tableau[i][noArtificials + noSlags + noVariables] < 0)
                    for (int j = 0; j < noArtificials + noSlags + noColumns; j++)
                        tableau[i][j] = -tableau[i][j];

                tableau[i][i] = 1;
            }

            int count = noConstraints;
            for (int j = 0; j < noVariables; j++)
            {
                if (boundaries[j].y != double.MaxValue)
                {
                    tableau[count][noArtificials + noSlags + j] = 1;
                    tableau[count][noArtificials + noSlags + noVariables] = boundaries[j].y - boundaries[j].x;
                    tableau[count][count] = 1;
                    tableau[count][count + noConstraints + countUpperLimit] = 1;
                    count++;
                }
            }

            s = reader.ReadLine();
            line = s.Split(',');
            c = 0;
            for (int j = 0; j < noVariables; j++)
            {
                zFunction[noArtificials + noSlags + j] = -double.Parse(line[j], CultureInfo.InvariantCulture);
            }
            zFunction[noArtificials + noSlags + noVariables] = -double.Parse(line[noVariables], CultureInfo.InvariantCulture);
            
            for (int i = 0; i < noRows; i++)
            {
                for (int j = noArtificials; j < noArtificials + noSlags + noColumns; j++)
                {
                    wFunction[j] -= tableau[i][j];
                }
            }

            reader.Close();

            PrintTableau();
            
            bool found = false;

            if (PhaseOne())
            {
                PrintTableau();
                found = PhaseTwo();
                // if found: zFunc[last] is max of org. zFunc
            }
                Console.WriteLine(found);
            PrintTableau();


                this.Close();
        }

        private bool PhaseOne()
        {
            int incommingBasis;
            int pivotRow;

            count = 0;
            while (true)
            {

                incommingBasis = FindIncommingBasis(0, wFunction);
                if (incommingBasis < 0)
                {
                    if (Math.Abs(wFunction[rightSideIndex]) < 0.00000001)
                        return true;
                    else
                        return false;
                }

                pivotRow = FindPivotRow(incommingBasis);
                if (pivotRow < 0)
                    return false;

                Pivoting(incommingBasis, pivotRow, 0);
                count++;
            }
        }

        private bool PhaseTwo()
        {
            int incommingBasis;
            int pivotRow;
            count = 0;
            while (true)
            {
                incommingBasis = FindIncommingBasis(noArtificials, zFunction);
                if (incommingBasis < 0)
                    return true;

                pivotRow = FindPivotRow(incommingBasis);
                if (pivotRow < 0)
                    return false;

                Pivoting(incommingBasis, pivotRow, noArtificials);
                count++;
            }
        }

        private int FindIncommingBasis(int startIndex, double[] objectFunction)
        {
            int incommingBasis = -1;

            for (int i = startIndex; i < rightSideIndex; i++)
            {
                if (objectFunction[i] < 0)
                {
                    incommingBasis = i;
                    break;
                }
            }
            return incommingBasis;
        }

        private int FindPivotRow(int incommingBasis)
        {
            int pivotRow = -1;

            double ratio;

            double min = double.MaxValue;

            for (int i = 0; i < noRows; i++)
            {
                if (tableau[i][incommingBasis] > 0)
                {
                    ratio = tableau[i][rightSideIndex] / tableau[i][incommingBasis];
                    if (ratio < min)
                    {
                        min = ratio;
                        pivotRow = i;
                    }
                }
            }

            return pivotRow;

        }

        private void Pivoting(int incomingBasis, int pivotRow, int startIndex)
        {
            double temp = tableau[pivotRow][incomingBasis];
            for (int i = startIndex; i <= rightSideIndex; i++)
            {
                tableau[pivotRow][i] /= temp;
            }

            for (int j = 0; j < noRows; j++)
            {
                if (j != pivotRow)
                {
                    temp = tableau[j][incomingBasis];
                    for (int i = 0; i <= rightSideIndex; i++)
                    {
                        {
                            tableau[j][i] -= tableau[pivotRow][i] * temp;
                            if (Math.Abs(tableau[j][i]) < epsilon)
                                tableau[j][i] = 0;
                        }
                    }
                }
            }

            temp = zFunction[incomingBasis];
            for (int i = startIndex; i <= rightSideIndex; i++)
            {
                zFunction[i] -= tableau[pivotRow][i] * temp;
                if (Math.Abs(zFunction[i]) < epsilon)
                    zFunction[i] = 0;
            }

            if (startIndex == 0)
            {
                temp = wFunction[incomingBasis];
                for (int i = 0; i <= rightSideIndex; i++)
                {
                    wFunction[i] -= tableau[pivotRow][i] * temp;
                    if (Math.Abs(wFunction[i]) < epsilon)
                        wFunction[i] = 0;
                }
            }
           // PrintTableau();
        }

        private void PrintTableau()
        {
            StreamWriter writer = new StreamWriter(@"C:\Users\JHee\Documents\MyData\Tab.txt", false);

            for (int i = 0; i < tableau.Length; i++)
            {
                for (int j = 0; j < tableau[0].Length; j++)
                    writer.Write(string.Format("{0,10:f2}", tableau[i][j]));

                writer.WriteLine();
            }
            writer.WriteLine();
            writer.Close();

        }
    }
}
