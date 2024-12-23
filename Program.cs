using System.Text;
using ConvertXyzForGaussian;

class Program
{
    static void Main(string[] args)
    {
        try
        {
            Console.Write("Informe a pasta de origem dos arquivos .xyz: ");
            string inputFolder = Console.ReadLine();

            Console.Write("Informe a pasta de destino para os arquivos .com: ");
            string outputFolder = Console.ReadLine();

            Console.Write("Informe o diretório no Ubuntu para mover os arquivos .com: ");
            string ubuntuPath = Console.ReadLine();

            if (!Directory.Exists(inputFolder))
            {
                Console.WriteLine("A pasta de origem não existe.");
                return;
            }

            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            if (!Directory.Exists(ubuntuPath))
            {
                Console.WriteLine("O diretório do Ubuntu não existe.");
                return;
            }

            string[] xyzFiles = Directory.GetFiles(inputFolder, "*.xyz");

            if (xyzFiles.Length == 0)
            {
                Console.WriteLine("Nenhum arquivo .xyz encontrado na pasta de origem.");
                Console.ReadKey();
                return;
            }

            foreach (var file in xyzFiles)
            {
                var lines = File.ReadAllLines(file);
                var atomos = new List<Molecula.Atomo>();

                for (int i = 2; i < lines.Length; i++) // Ignorar as duas primeiras linhas
                {
                    var parts = lines[i].Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    if (parts.Length == 4)
                    {
                        atomos.Add(new Molecula.Atomo
                        {
                            Elemento = parts[0],
                            X = double.Parse(parts[1]),
                            Y = double.Parse(parts[2]),
                            Z = double.Parse(parts[3])
                        });
                    }
                }

                if (atomos.Count >= 4) // Ao menos 4 átomos para calcular diedro
                {
                    var molecula = new Molecula();

                    StringBuilder gaussianFile = new StringBuilder();
                    gaussianFile.AppendLine("%Mem=8GB");
                    gaussianFile.AppendLine("%NProcShared=8");
                    gaussianFile.AppendLine("#P HF/6-31G* Opt");
                    gaussianFile.AppendLine("");
                    gaussianFile.AppendLine(Path.GetFileNameWithoutExtension(file));
                    gaussianFile.AppendLine("0 1");

                    for (int i = 0; i < atomos.Count; i++)
                    {
                        var atomo = atomos[i];
                        gaussianFile.Append($"{atomo.Elemento} ");

                        if (i >= 1)
                        {
                            double distancia = molecula.CalcularDistancia(atomos, i, i - 1);
                            gaussianFile.Append($"B{i} {distancia:F6} ");
                        }

                        if (i >= 2)
                        {
                            double angulo = molecula.CalcularAngulo(atomos, i, i - 1, i - 2);
                            gaussianFile.Append($"A{i} {angulo:F6} ");
                        }

                        if (i >= 3)
                        {
                            double diedro = molecula.CalcularDiedro(atomos, i, i - 1, i - 2, i - 3);
                            gaussianFile.Append($"D{i} {diedro:F6} ");
                        }

                        gaussianFile.AppendLine();
                    }

                    gaussianFile.AppendLine("");

                    string outputFile = Path.Combine(outputFolder, Path.GetFileNameWithoutExtension(file) + ".com");
                    File.WriteAllText(outputFile, gaussianFile.ToString());
                }
            }

            MoveComFilesToUbuntu(outputFolder, ubuntuPath);
            Console.WriteLine("Conversão concluída. Pressione qualquer tecla para sair.");
            Console.ReadKey();
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Erro: {ex.Message}");
        }
    }

    static void MoveComFilesToUbuntu(string localSavePath, string ubuntuPath)
    {
        string[] files = Directory.GetFiles(localSavePath, "*.com");
        foreach (string file in files)
        {
            string fileName = Path.GetFileName(file);
            string destFile = Path.Combine(ubuntuPath, fileName);
            try
            {
                File.Move(file, destFile);
                Console.WriteLine($"Movido {fileName} para {ubuntuPath}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Erro ao mover {fileName}: {ex.Message}");
            }
        }
    }
}

namespace ConvertXyzForGaussian
{
    public class Molecula
    {
        public double CalcularDistancia(List<Atomo> atomos, int i, int j)
        {
            var atomoI = atomos[i];
            var atomoJ = atomos[j];
            return Math.Sqrt(
                Math.Pow(atomoJ.X - atomoI.X, 2) +
                Math.Pow(atomoJ.Y - atomoI.Y, 2) +
                Math.Pow(atomoJ.Z - atomoI.Z, 2)
            );
        }

        public double CalcularAngulo(List<Atomo> atomos, int i, int j, int k)
        {
            var vecIJ = Vetor(atomos[i], atomos[j]);
            var vecJK = Vetor(atomos[j], atomos[k]);
            return CalcularAnguloEntreVetores(vecIJ, vecJK);
        }

        public double CalcularDiedro(List<Atomo> atomos, int i, int j, int k, int l)
        {
            var vecIJ = Vetor(atomos[i], atomos[j]);
            var vecJK = Vetor(atomos[j], atomos[k]);
            var vecKL = Vetor(atomos[k], atomos[l]);
            return CalcularAnguloEntreVetores(CrossProduct(vecIJ, vecJK), CrossProduct(vecJK, vecKL));
        }

        private double[] Vetor(Atomo a, Atomo b) => new[] { b.X - a.X, b.Y - a.Y, b.Z - a.Z };

        private double[] CrossProduct(double[] a, double[] b) => new[]
        {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        };

        private double CalcularAnguloEntreVetores(double[] a, double[] b)
        {
            double produtoEscalar = a.Zip(b, (x, y) => x * y).Sum();
            double magnitudeA = Math.Sqrt(a.Sum(x => x * x));
            double magnitudeB = Math.Sqrt(b.Sum(x => x * x));
            return Math.Acos(produtoEscalar / (magnitudeA * magnitudeB)) * (180 / Math.PI);
        }

        public class Atomo
        {
            public string Elemento { get; set; }
            public double X { get; set; }
            public double Y { get; set; }
            public double Z { get; set; }
        }
    }
}




// using System.Text;
// using ConvertXyzForGaussian;
//
// class Program
// {
//     static void Main(string[] args)
//     {
//         try
//         {
//             Console.Write("Informe a pasta de origem dos arquivos .xyz: ");
//             string inputFolder = Console.ReadLine();
//
//             Console.Write("Informe a pasta de destino para os arquivos .com: ");
//             string outputFolder = Console.ReadLine();
//
//             Console.Write("Informe o diretório no Ubuntu para mover os arquivos .com: ");
//             string ubuntuPath = Console.ReadLine();
//
//             if (!Directory.Exists(inputFolder))
//             {
//                 Console.WriteLine("A pasta de origem não existe.");
//                 return;
//             }
//
//             if (!Directory.Exists(outputFolder))
//             {
//                 Directory.CreateDirectory(outputFolder);
//             }
//
//             if (!Directory.Exists(ubuntuPath))
//             {
//                 Console.WriteLine("O diretório do Ubuntu não existe.");
//                 return;
//             }
//
//             string[] xyzFiles = Directory.GetFiles(inputFolder, "*.xyz");
//
//             if (xyzFiles.Length == 0)
//             {
//                 Console.WriteLine("Nenhum arquivo .xyz encontrado na pasta de origem.");
//                 Console.WriteLine("Precione qualquer tecla para sair.");
//                 Console.ReadKey();
//                 Environment.Exit(0);
//                 return;
//             }
//             
//             //Infome as dimençoes 
//             int distanciA = 0;
//             int distanciB = 1;
//             
//             int anguloA = 1;
//             int anguloB = 0;
//             int anguloC = 2;
//             
//             int diedroA = 1;
//             int diedroB = 0;
//             int diedroC = 2;
//             int diedroD = 3;
//             
//
//             foreach (var file in xyzFiles)
//             {
//                 var lines = File.ReadAllLines(file);
//                 var atomos = new List<Molecula.Atomo>();
//
//                 // Lê o arquivo a partir da terceira linha
//                 for (int i = 2; i < lines.Length; i++)
//                 {
//                     var calc = lines[i].Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
//                     if (calc.Length == 4) // Verifica se a linha tem os elementos esperados
//                     {
//                         var atomo = new Molecula.Atomo
//                         {
//                             Elemento = calc[0],
//                             X = double.Parse(calc[1]),
//                             Y = double.Parse(calc[2]),
//                             Z = double.Parse(calc[3])
//                         };
//
//                         //Função calcula o Gaussion
//                         var funcao = new Molecula();
//                         atomo.X = funcao.CalcularDistancia(atomo);
//                         atomo.Y = funcao.CalcularAngulo(atomo);
//                         atomo.Z = funcao.CalcularDiedro(atomo);
//
//                         atomos.Add(atomo);
//                     }
//                 }
//
//                 if (atomos.Any())
//                 {
//                     // Aqui é gerado o conteúdo no formato .com
//                     StringBuilder gaussianFile = new StringBuilder(outputFolder);
//                     gaussianFile.AppendLine("%Mem=8GB");
//                     gaussianFile.AppendLine("%NProcShared=8");
//                     gaussianFile.AppendLine("#P HF/6-31G* Opt");
//                     gaussianFile.AppendLine("");
//                     gaussianFile.AppendLine("");
//                     gaussianFile.AppendLine("0 1"); // Multiplicidade e carga
//                     foreach (var line in lines)
//                     {
//                         if (line.Trim().Length > 0 && char.IsLetter(line[0])) // Evita linhas em branco
//                         {
//                             gaussianFile.AppendLine(line);
//                         }
//                     }
//
//                     gaussianFile.AppendLine("");
//                 }
//                 
//             }
//
//             //Mover para a pasta do ubuntu.
//             MoveComFilesToUbuntu(outputFolder, ubuntuPath);
//             Console.WriteLine("Conversão Concluído.");
//             Console.WriteLine("Preciona qualquer tecla p/ finalizar o programa...");
//             Console.ReadLine();
//         }
//         catch (Exception ex)
//         {
//             Console.WriteLine($"Erro: {ex.Message}");
//         }
//     }
//
//     static void MoveComFilesToUbuntu(string localSavePath, string ubuntuPath)
//     {
//         string[] files = Directory.GetFiles(localSavePath, "*.com");
//         foreach (string file in files)
//         {
//             string fileName = Path.GetFileName(file);
//             string destFile = Path.Combine(ubuntuPath, fileName);
//             try
//             {
//                 File.Move(file, destFile);
//                 Console.WriteLine($"Movido {fileName} para {ubuntuPath}");
//             }
//             catch (Exception ex)
//             {
//                 Console.WriteLine($"Erro ao mover {fileName}: {ex.Message}");
//             }
//         }
//     }
// }