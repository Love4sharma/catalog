import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import org.json.JSONObject;

public class ShamirSecretMatrix {

    
    static class Point {
        int x;
        BigInteger y;

        Point(int x, BigInteger y) {
            this.x = x;
            this.y = y;
        }
    }

    
    public static BigInteger decodeValue(String value, int base) {
        return new BigInteger(value, base);
    }

    /
    public static BigInteger[] gaussianElimination(double[][] A, double[] Y) {
        int n = A.length;

        for (int i = 0; i < n; i++) {
            double maxEl = Math.abs(A[i][i]);
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(A[k][i]) > maxEl) {
                    maxEl = Math.abs(A[k][i]);
                    maxRow = k;
                }
            }

            for (int k = i; k < n; k++) {
                double tmp = A[maxRow][k];
                A[maxRow][k] = A[i][k];
                A[i][k] = tmp;
            }
            double tmp = Y[maxRow];
            Y[maxRow] = Y[i];
            Y[i] = tmp;

            for (int k = i + 1; k < n; k++) {
                double c = -A[k][i] / A[i][i];
                for (int j = i; j < n; j++) {
                    if (i == j) {
                        A[k][j] = 0;
                    } else {
                        A[k][j] += c * A[i][j];
                    }
                }
                Y[k] += c * Y[i];
            }
        }

        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            x[i] = Y[i] / A[i][i];
            for (int k = i - 1; k >= 0; k--) {
                Y[k] -= A[k][i] * x[i];
            }
        }

        BigInteger[] result = new BigInteger[n];
        for (int i = 0; i < n; i++) {
            result[i] = BigInteger.valueOf(Math.round(x[i]));
        }

        return result;
    }

    
    public static double[][] constructVandermondeMatrix(List<Point> points, int degree) {
        int k = points.size();
        double[][] A = new double[k][degree + 1];
        for (int i = 0; i < k; i++) {
            int x = points.get(i).x;
            for (int j = 0; j <= degree; j++) {
                A[i][degree - j] = Math.pow(x, j);
            }
        }
        return A;
    }

   
    public static double[] extractYValues(List<Point> points) {
        double[] Y = new double[points.size()];
        for (int i = 0; i < points.size(); i++) {
            Y[i] = points.get(i).y.doubleValue();
        }
        return Y;
    }

    
    public static BigInteger findSecret(List<Point> points) {
        int degree = points.size() - 1;
        double[][] A = constructVandermondeMatrix(points, degree);
        double[] Y = extractYValues(points);

        BigInteger[] coefficients = gaussianElimination(A, Y);

        return coefficients[degree];
    }

   
    public static String readJsonFromFile(String filePath) throws IOException {
        return new String(Files.readAllBytes(Paths.get(filePath)));
    }

   
    public static List<Point> parseJsonData(String jsonData) {
        JSONObject json = new JSONObject(jsonData);
        List<Point> points = new ArrayList<>();
        JSONObject keys = json.getJSONObject("keys");
        int n = keys.getInt("n");
        int k = keys.getInt("k");

        for (String key : json.keySet()) {
            if (!key.equals("keys")) {
                JSONObject root = json.getJSONObject(key);
                int x = Integer.parseInt(key); 
                int base = root.getInt("base"); 
                String value = root.getString("value"); 
                BigInteger y = decodeValue(value, base); 
                points.add(new Point(x, y));
            }
        }
        return points;
    }

    public static void main(String[] args) {
        try {
            
            String filePath = "src/data.json";
            // System.out.println(filePath);

            
            String jsonData = readJsonFromFile(filePath);

            
            List<Point> points = parseJsonData(jsonData);

            
            BigInteger secret = findSecret(points);

            //  constant term 'c'
            System.out.println("The secret constant term (c) is: " + secret);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
