namespace CourseworkLibrary
{
    public class Results
    {
        public static void WriteToFile(double[,,] solution)
        {
            string filename = "Result.txt";
            try
            {
                using (StreamWriter writer = new StreamWriter(filename))
                {
                    int nx = solution.GetLength(0);
                    int ny = solution.GetLength(1);
                    int nz = solution.GetLength(2);

                    for (int i = 0; i < nx; i++)
                    {
                        for (int j = 0; j < ny; j++)
                        {
                            for (int k = 0; k < nz; k++)
                            {
                                writer.WriteLine($"{i} {j} {k} {solution[i, j, k]}"); // Сохраняем координаты и температуру
                            }
                        }
                    }
                }
                Console.Clear();
                Console.WriteLine($"Решение записано в {filename}");
            }
            catch (Exception e)
            {
                Console.Clear();
                Console.WriteLine($"Error writing to file: {e.Message}");
            }
        }

        public static void PrintSlise(double[,,] solution)
        {
            while (true)
            {
                Console.Clear();
                Console.WriteLine("\n1) Выберите ось:\nX\nY\nZ\n\nQ - Назад");
                string axis = Console.ReadLine();
                if (axis == "q" || axis == "Q") break;
                int size = axis switch
                {
                    "x" or "X" => solution.GetLength(0),
                    "y" or "Y" => solution.GetLength(1),
                    "z" or "Z" => solution.GetLength(2),
                    _ => throw new ArgumentException("Недопустимая ось. Используйте 'x', 'y' или 'z'")
                };
                Console.Clear();
                Console.WriteLine($"\nПослойные срезы по оси {axis.ToUpper()}:");
                Thread.Sleep(1000);
                Console.WriteLine("(Для навигации используйте стрелки, для выхода - Esc)\n");
                Thread.Sleep(2000);

                int currentSlice = 0;
                ConsoleKey key;

                while (true)
                {
                    Console.Clear();

                    PrintSingleSlice(solution, axis, currentSlice);

                    Console.WriteLine($"\nСрез {axis.ToUpper()} = {currentSlice} (Всего: {size - 1})");
                    Console.WriteLine("<- ->  - листать слои | | Esc - выход");

                    key = Console.ReadKey(true).Key;
                    if (key == ConsoleKey.Escape)
                        break;

                    currentSlice = key switch
                    {
                        ConsoleKey.LeftArrow => currentSlice == 0 ? currentSlice : currentSlice - 1,
                        ConsoleKey.RightArrow => currentSlice == size - 1 ? currentSlice : currentSlice + 1
                    };

                }
            }
        }

        private static void PrintSingleSlice(double[,,] solution, string axis, int slice)
        {
            switch (axis.ToLower())
            {
                case "x":
                    PrintSliceYZ(solution, slice);
                    break;
                case "y":
                    PrintSliceXZ(solution, slice);
                    break;
                case "z":
                    PrintSliceXY(solution, slice);
                    break;
            }
        }

        private static void PrintSliceXY(double[,,] solution, int z)
        {
            Console.WriteLine($"┌{new string('─', solution.GetLength(1) * 10 + 1)}┐");
            Console.Write("│   X→ ");
            for (int x = 0; x < solution.GetLength(0); x++) Console.Write($"{x,8}");
            Console.WriteLine(" │");
            Console.WriteLine($"├{new string('─', solution.GetLength(1) * 10 + 1)}┤");

            for (int y = 0; y < solution.GetLength(1); y++)
            {
                Console.Write($"│ Y{y,2} │");
                for (int x = 0; x < solution.GetLength(0); x++)
                {
                    Console.Write($"{solution[x, y, z],8:F2}");
                }
                Console.WriteLine(" │");
            }
            Console.WriteLine($"└{new string('─', solution.GetLength(1) * 10 + 1)}┘");
        }

        private static void PrintSliceXZ(double[,,] solution, int y)
        {
            Console.WriteLine($"┌{new string('─', solution.GetLength(0) * 10 + 1)}┐");
            Console.Write("│   X→ ");
            for (int x = 0; x < solution.GetLength(0); x++) Console.Write($"{x,8}");
            Console.WriteLine(" │");
            Console.WriteLine($"├{new string('─', solution.GetLength(0) * 10 + 1)}┤");

            for (int z = 0; z < solution.GetLength(2); z++)
            {
                Console.Write($"│ Z{z,2} │");
                for (int x = 0; x < solution.GetLength(0); x++)
                {
                    Console.Write($"{solution[x, y, z],8:F2}");
                }
                Console.WriteLine(" │");
            }
            Console.WriteLine($"└{new string('─', solution.GetLength(0) * 10 + 1)}┘");
        }

        private static void PrintSliceYZ(double[,,] solution, int x)
        {
            Console.WriteLine($"┌{new string('─', solution.GetLength(1) * 10 + 1)}┐");
            Console.Write("│   Y→ ");
            for (int y = 0; y < solution.GetLength(1); y++) Console.Write($"{y,8}");
            Console.WriteLine(" │");
            Console.WriteLine($"├{new string('─', solution.GetLength(1) * 10 + 1)}┤");

            for (int z = 0; z < solution.GetLength(2); z++)
            {
                Console.Write($"│ Z{z,2} │");
                for (int y = 0; y < solution.GetLength(1); y++)
                {
                    Console.Write($"{solution[x, y, z],8:F2}");
                }
                Console.WriteLine(" │");
            }
            Console.WriteLine($"└{new string('─', solution.GetLength(1) * 10 + 1)}┘");
        }
    }
}
