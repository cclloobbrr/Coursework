using MathNet.Numerics.LinearAlgebra;

namespace PracticeLibrary
{
    public class Create
    {
        public static Matrix<double> TridiagonalMatrixX(int nx, double rx)
        {
            var matrix = Matrix<double>.Build.Dense(nx, nx, 0.0);

            for (int i = 1; i < nx - 1; i++)
            {
                matrix[i, i - 1] = -rx;
                matrix[i, i] = 1 + 2 * rx;
                matrix[i, i + 1] = -rx;
            }

            return matrix;
        }

        public static Matrix<double> TridiagonalMatrixY(int ny, double ry)
        {
            var matrix = Matrix<double>.Build.Dense(ny, ny, 0.0);

            for (int i = 1; i < ny - 1; i++)
            {
                matrix[i, i - 1] = -ry;
                matrix[i, i] = 1 + 2 * ry;
                matrix[i, i + 1] = -ry;
            }

            return matrix;
        }

        public static Matrix<double> TridiagonalMatrixZ(int nz, double rz)
        {
            var matrix = Matrix<double>.Build.Dense(nz, nz, 0.0);

            for (int i = 1; i < nz - 1; i++)
            {
                matrix[i, i - 1] = -rz;
                matrix[i, i] = 1 + 2 * rz;
                matrix[i, i + 1] = -rz;
            }

            return matrix;
        }

        public static Vector<double> RightHandSideX(
            double[,,] solution,
            double alpha,
            double dt,
            double rx,
            int nx,
            int j,
            int k,
            double dx,
            double dy,
            double dz,
            double t)
        {
            Vector<double> vector = Vector<double>.Build.Dense(nx, 0.0);

            for (int i = 1; i < nx - 1; i++)
            {
                vector[i] = solution[i, j, k] +
                          rx * (solution[i - 1, j, k] - 2 * solution[i, j, k] + solution[i + 1, j, k]);
            }

            return vector;
        }

        public static Vector<double> RightHandSideY(
            double[,,] solution,
            double alpha,
            double dt,
            double ry,
            int ny,
            int i,
            int k,
            double dx,
            double dy,
            double dz,
            double t)
        {
            Vector<double> vector = Vector<double>.Build.Dense(ny, 0.0);

            for (int j = 1; j < ny - 1; j++)
            {
                vector[j] = solution[i, j, k] +
                          ry * (solution[i, j - 1, k] - 2 * solution[i, j, k] + solution[i, j + 1, k]);
            }

            return vector;
        }

        public static Vector<double> RightHandSideZ(
            double[,,] solution,
            double alpha,
            double dt,
            double rz,
            int nz,
            int i,
            int j,
            double dx,
            double dy,
            double dz,
            double t)
        {
            Vector<double> vector = Vector<double>.Build.Dense(nz, 0.0);

            for (int k = 1; k < nz - 1; k++)
            {
                vector[k] = solution[i, j, k] +
                          rz * (solution[i, j, k - 1] - 2 * solution[i, j, k] + solution[i, j, k + 1]);
            }

            return vector;
        }
    }
}
