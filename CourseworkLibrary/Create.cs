using MathNet.Numerics.LinearAlgebra;

namespace CourseworkLibrary
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

        public static Vector<double> RightHandSideX(double[,,] solution, double alpha, double dt, double rx, int nx, int j, int k, double dx, double dy, double dz, List<(double x, double y, double z, double Q)> pointSources, double t)
        {
            Vector<double> vector = Vector<double>.Build.Dense(nx, 0.0);
            double y = j * dy;
            double z = k * dz;
            for (int i = 1; i < nx - 1; i++)
            {
                double x = i * dx;
                vector[i] = solution[i, j, k] + rx * (solution[i - 1, j, k] - 2 * solution[i, j, k] + solution[i + 1, j, k]);

                if (pointSources != null)
                {
                    foreach (var source in pointSources)
                    {
                        double distance = Math.Sqrt(Math.Pow(x - source.x, 2) + Math.Pow(y - source.y, 2) + Math.Pow(z - source.z, 2));
                        if (distance <= dx / 2)
                        {
                            vector[i] += alpha * dt * source.Q / (dx * dy * dz);
                        }
                    }
                }
            }
            return vector;
        }

        public static Vector<double> RightHandSideY(double[,,] solution, double alpha, double dt, double ry, int ny, int i, int k, double dx, double dy, double dz, List<(double x, double y, double z, double Q)> pointSources, double t)
        {
            Vector<double> vector = Vector<double>.Build.Dense(ny, 0.0);

            double x = i * dx;
            double z = k * dz;

            for (int j = 1; j < ny - 1; j++)
            {
                double y = j * dy;

                vector[j] = solution[i, j, k] + ry * (solution[i, j - 1, k] - 2 * solution[i, j, k] + solution[i, j + 1, k]);

                if (pointSources != null)
                {
                    foreach (var source in pointSources)
                    {
                        double distance = Math.Sqrt(Math.Pow(x - source.x, 2) + Math.Pow(y - source.y, 2) + Math.Pow(z - source.z, 2));
                        if (distance <= dy / 2)
                        {
                            vector[j] += alpha * dt * source.Q / (dx * dy * dz);
                        }
                    }
                }
            }
            return vector;
        }

        public static Vector<double> RightHandSideZ(double[,,] solution, double alpha, double dt, double rz, int nz, int i, int j, double dx, double dy, double dz, List<(double x, double y, double z, double Q)> pointSources, double t)
        {
            Vector<double> vector = Vector<double>.Build.Dense(nz, 0.0);
            double x = i * dx;
            double y = j * dy;

            for (int k = 1; k < nz - 1; k++)
            {
                double z = k * dz;
                vector[k] = solution[i, j, k] + rz * (solution[i, j, k - 1] - 2 * solution[i, j, k] + solution[i, j, k + 1]);

                if (pointSources != null)
                {
                    foreach (var source in pointSources)
                    {
                        double distance = Math.Sqrt(Math.Pow(x - source.x, 2) + Math.Pow(y - source.y, 2) + Math.Pow(z - source.z, 2));

                        if (distance <= dz / 2)
                        {
                            vector[k] += alpha * dt * source.Q / (dx * dy * dz);
                        }
                    }
                }
            }
            return vector;
        }
    }
}
