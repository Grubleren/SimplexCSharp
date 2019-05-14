using System;
using System.Windows.Forms;
using System.IO;
using System.Drawing;
using System.Globalization;
using System.Collections.Generic;

namespace SimplexCSharp
{
    public partial class Form1 : Form
    {
        double epsilon = 0.0000000000001;
        new struct Bounds
        {
            public double lower;
            public double upper;
        }

        int count;

        int nConstraints;
        int nVariables;
        int nRows;
        int nSlags;
        int nArtificials;
        int slagIndex;
        int varibleIndex;
        int righthandSideIndex;

        double[][] tableau;
        double[] zFunction;
        double[] wFunction;
        int[] basis;
        double[] solution;

        bool[] allowedRows;
        bool[] allowedColumns;

        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            string inputDataPath = @"C:..\..\TestData5.txt";
            StreamReader reader = new StreamReader(inputDataPath);

            string s = SkipComments(reader);

            ProblemSize(s);

            Bounds[] boundaries = new Bounds[nVariables];

            GetBoundaries(reader, boundaries);

            AllocateTableauAndObjectfunction();

            AddConstraints(reader, boundaries);

            AddUpperLimitConstraints(reader, boundaries);

            AddObjectFunction(reader);

            allowedRows = new bool[nRows];
            allowedColumns = new bool[righthandSideIndex + 1];

            RemoveUnbounded(boundaries);

            AddArtificials();

            AddArtificialObjectFunction(reader);

            reader.Close();

            PrintTableau();

            bool found = false;

            if (PhaseOne())
            {
                PrintTableau();
                found = PhaseTwo();
                // if found: zFunc[last] is max of org. zFunc
            }

            if (found)
            {
                for (int i = 0; i < righthandSideIndex; i++)
                {
                    if (allowedColumns[i])
                    {
                        int row = InBasis(i);
                        double solu;
                        if (row >= 0)
                            solu = tableau[row][righthandSideIndex];
                        else
                            solu = 0;
                        solution[i] = solu;
                    }
                }

                for (int i = varibleIndex; i < righthandSideIndex; i++)
                {
                    if (!allowedColumns[i])
                    {
                        int row = InBasis(i);
                        double solu = tableau[row][righthandSideIndex];
                        for (int j = 0; j < righthandSideIndex; j++)
                        {
                            if (allowedColumns[j])
                                solu -= tableau[row][j] * solution[j];

                        }
                        solution[i] = solu;
                    }
                }

                for (int i = varibleIndex; i < righthandSideIndex; i++)
                {
                    if (allowedColumns[i])
                    {
                        solution[i] += boundaries[i - varibleIndex].lower;
                    }
                }


                Console.WriteLine("Artificials:");
                for (int i = 0; i < slagIndex; i++)
                    Console.WriteLine("{0:0.######}", solution[i]);
                Console.WriteLine("Slags:");
                for (int i = slagIndex; i < varibleIndex; i++)
                    Console.WriteLine("{0:0.######}", solution[i]);
                Console.WriteLine("Variables");
                for (int i = varibleIndex; i < righthandSideIndex; i++)
                    Console.WriteLine("{0:0.######}", solution[i]);
            }
            else
                Console.WriteLine("No solution");
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
                    if (Math.Abs(wFunction[righthandSideIndex]) < 0.00000001)
                        return true;
                    else
                        return false;
                }

                pivotRow = FindPivotRow(incommingBasis);
                if (pivotRow < 0)
                    return false;

                basis[pivotRow] = incommingBasis;

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
                incommingBasis = FindIncommingBasis(nArtificials, zFunction);
                if (incommingBasis < 0)
                    return true;

                pivotRow = FindPivotRow(incommingBasis);
                if (pivotRow < 0)
                    return false;

                basis[pivotRow] = incommingBasis;

                Pivoting(incommingBasis, pivotRow, nArtificials);
                count++;
            }
        }

        private int FindIncommingBasis(int startIndex, double[] objectFunction)
        {
            int incommingBasis = -1;

            for (int i = startIndex; i < righthandSideIndex; i++)
            {
                if (allowedColumns[i])
                {
                    if (objectFunction[i] < 0)
                    {
                        incommingBasis = i;
                        break;
                    }
                }
            }
            return incommingBasis;
        }

        private int FindPivotRow(int incommingBasis)
        {
            int pivotRow = -1;

            double ratio;

            double min = double.MaxValue;

            for (int i = 0; i < nRows; i++)
            {
                if (allowedRows[i])
                {
                    if (tableau[i][incommingBasis] > 0)
                    {
                        ratio = tableau[i][righthandSideIndex] / tableau[i][incommingBasis];
                        if (ratio < min)
                        {
                            min = ratio;
                            pivotRow = i;
                        }
                    }
                }
            }

            return pivotRow;

        }

        private void Pivoting(int incomingBasis, int pivotRow, int startIndex)
        {
            double temp = tableau[pivotRow][incomingBasis];
            for (int i = startIndex; i <= righthandSideIndex; i++)
            {
                if (allowedColumns[i])
                    tableau[pivotRow][i] /= temp;
            }

            for (int j = 0; j < nRows; j++)
            {
                if (allowedRows[j])
                {
                    if (j != pivotRow)
                    {
                        temp = tableau[j][incomingBasis];
                        for (int i = 0; i <= righthandSideIndex; i++)
                        {
                            {
                                tableau[j][i] -= tableau[pivotRow][i] * temp;
                                if (Math.Abs(tableau[j][i]) < epsilon)
                                    tableau[j][i] = 0;
                            }
                        }
                    }
                }
            }

            temp = zFunction[incomingBasis];
            for (int i = startIndex; i <= righthandSideIndex; i++)
            {
                if (allowedColumns[i])
                {
                    zFunction[i] -= tableau[pivotRow][i] * temp;
                    if (Math.Abs(zFunction[i]) < epsilon)
                        zFunction[i] = 0;
                }
            }

            if (startIndex == 0)
            {
                temp = wFunction[incomingBasis];
                for (int i = 0; i <= righthandSideIndex; i++)
                {
                    if (allowedColumns[i])
                    {
                        wFunction[i] -= tableau[pivotRow][i] * temp;
                        if (Math.Abs(wFunction[i]) < epsilon)
                            wFunction[i] = 0;
                    }
                }
            }
            // PrintTableau();
        }

        private void PrintTableau()
        {
            StreamWriter writer = new StreamWriter(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile) + @"\Tableau.txt", false);

            for (int i = 0; i < tableau.Length; i++)
            {
                for (int j = 0; j < tableau[0].Length; j++)
                    writer.Write(string.Format("{0,10:f2}", tableau[i][j]));

                writer.WriteLine();
            }
            writer.WriteLine();
            writer.Close();

        }

        string SkipComments(StreamReader reader)
        {
            string s = reader.ReadLine();
            while (s[0] == '#')
                s = reader.ReadLine();
            return s;
        }

        void ProblemSize(string s)
        {
            string[] size = s.Split(',');
            nConstraints = int.Parse(size[0]);
            nVariables = int.Parse(size[1]);
        }

        void GetBoundaries(StreamReader reader, Bounds[] boundaries)
        {
            string s = reader.ReadLine();
            string[] limits = s.Split(',');

            int countUpperLimit = 0;
            for (int i = 0; i < nVariables; i++)
            {
                boundaries[i].lower = !limits[2 * i].Contains("inf") ? double.Parse(limits[2 * i], CultureInfo.InvariantCulture) : double.MinValue;
                boundaries[i].upper = !limits[2 * i + 1].Contains("inf") ? double.Parse(limits[2 * i + 1], CultureInfo.InvariantCulture) : double.MaxValue;
                if (boundaries[i].upper != double.MaxValue)
                    countUpperLimit++;
            }

            nRows = nConstraints + countUpperLimit;
            nSlags = nRows;
            nArtificials = nRows;
            slagIndex = nArtificials;
            varibleIndex = slagIndex + nSlags;
            righthandSideIndex = nArtificials + nSlags + nVariables;
        }

        void AllocateTableauAndObjectfunction()
        {
            tableau = new double[nRows][];
            for (int i = 0; i < nRows; i++)
                tableau[i] = new double[righthandSideIndex + 1];

            zFunction = new double[righthandSideIndex + 1];
            wFunction = new double[righthandSideIndex + 1];
            solution = new double[righthandSideIndex];
            basis = new int[righthandSideIndex];

        }

        void AddConstraints(StreamReader reader, Bounds[] boundaries)
        {
            for (int i = 0; i < nConstraints; i++)
            {
                string s = reader.ReadLine();
                string[] line = s.Split(',');
                double lowerLimitSum = 0;
                for (int j = 0; j < nVariables; j++)
                {
                    double d = double.Parse(line[j], CultureInfo.InvariantCulture);
                    tableau[i][varibleIndex + j] = d;
                    lowerLimitSum += boundaries[j].lower != double.MinValue ? boundaries[j].lower * d : 0;
                }
                tableau[i][righthandSideIndex] = double.Parse(line[nVariables + 1], CultureInfo.InvariantCulture) - lowerLimitSum;

                if (line[nVariables].Contains("<"))
                    tableau[i][nArtificials + i] = 1;
                if (line[nVariables].Contains(">"))
                    tableau[i][nArtificials + i] = -1;
                if (tableau[i][righthandSideIndex] < 0)
                    for (int j = 0; j < righthandSideIndex + 1; j++)
                        tableau[i][j] = -tableau[i][j];
            }
        }

        void AddUpperLimitConstraints(StreamReader reader, Bounds[] boundaries)
        {
            int count = nConstraints;
            for (int j = 0; j < nVariables; j++)
            {
                if (boundaries[j].upper != double.MaxValue)
                {
                    tableau[count][varibleIndex + j] = 1;
                    tableau[count][righthandSideIndex] = boundaries[j].upper - (boundaries[j].lower == double.MinValue ? 0 :  boundaries[j].lower);
                    tableau[count][count + nArtificials] = 1;
                    count++;
                }
            }
        }

        void AddObjectFunction(StreamReader reader)
        {
            string s = reader.ReadLine();
            string[] line = s.Split(',');
            for (int j = 0; j < nVariables; j++)
            {
                zFunction[varibleIndex + j] = -double.Parse(line[j], CultureInfo.InvariantCulture);
            }
            zFunction[righthandSideIndex] = -double.Parse(line[nVariables], CultureInfo.InvariantCulture);

        }

        void RemoveUnbounded(Bounds[] boundaries)
        {
            List<int> unbounded = new List<int>();
            for (int i = 0; i < righthandSideIndex + 1; i++)
                allowedColumns[i] = true;

            for (int i = 0; i < nVariables; i++)
            {
                if (boundaries[i].lower == double.MinValue)
                {
                    unbounded.Add(varibleIndex + i);
                    allowedColumns[varibleIndex + i] = false;
                }
            }
            
            List<int> freeRows = new List<int>();
            for (int j = 0; j < nRows; j++)
            {
                freeRows.Add(j);
                allowedRows[j] = true;
            }
            for (int n = 0; n < unbounded.Count; n++)
            {
                int i = unbounded[n];
                double max = 0;
                int maxj = 0;
                foreach (int j in freeRows)
                {
                    if (Math.Abs(tableau[j][i]) > max)
                    {
                        max = Math.Abs(tableau[j][i]);
                        maxj = j;
                    }
                }
                freeRows.Remove(maxj);
                allowedRows[maxj] = false;
                basis[maxj] = i;
                double pivotElement = tableau[maxj][i];
                for (int k = slagIndex; k < righthandSideIndex + 1; k++)
                    tableau[maxj][k] /= pivotElement;

                for (int j = 0; j < nRows; j++)
                {
                    if (j != maxj)
                    {
                        pivotElement = tableau[j][i];
                        for (int k = slagIndex; k < righthandSideIndex + 1; k++)
                            tableau[j][k] -= tableau[maxj][k] * pivotElement;

                        if (tableau[j][righthandSideIndex] < 0)
                            for (int k = slagIndex; k < righthandSideIndex + 1; k++)
                                tableau[j][k] *= -1;
                    }
                }

                pivotElement = zFunction[i];
                for (int k = slagIndex; k < righthandSideIndex + 1; k++)
                    zFunction[k] -= tableau[maxj][k] * pivotElement;

            }
        }

        void AddArtificials()
        {
            for (int i = 0; i < nRows; i++)
                tableau[i][i] = 1;

            for (int i = 0; i < nRows; i++)
                if (allowedRows[i])
                    basis[i] = i;
        }

        void AddArtificialObjectFunction(StreamReader reader)
        {

            for (int i = 0; i < nRows; i++)
            {
                if (allowedRows[i])
                {
                    for (int j = slagIndex; j < righthandSideIndex + 1; j++)
                    {
                        if (allowedColumns[j])
                            wFunction[j] -= tableau[i][j];
                    }
                }
            }
        }

        int InBasis(int column)
        {
            int row = -1;
            for (int i = 0; i < nRows; i++)
                 if (basis[i] == column)
                {
                    row = i;
                    break;
                }
            return row;
        }
    }
}
