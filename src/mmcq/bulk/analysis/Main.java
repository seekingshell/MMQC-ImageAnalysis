package mmcq.bulk.analysis;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

import mmcq.bulk.analysis.MMCQ.VBox;
import mmcq.bulk.analysis.MMCQ.CMap;

public class Main {
    public static void main(String args[]) throws IOException {
        List<String> imageUrls = Files.readAllLines(Paths.get("urls.txt"));

        PrintWriter writer = new PrintWriter(new File("MMCQ_image_analysis.csv"));
        for(String imageUrl : imageUrls) {
            ArrayList<VBox> sortedColorPalette = processImage(imageUrl);
            int topOfList = sortedColorPalette.size() - 1;

            StringBuilder output = new StringBuilder();
            output.append(imageUrl + ",");
            output.append(formatRGBColor(sortedColorPalette.get(topOfList)) + ",");
            output.append(formatRGBColor(sortedColorPalette.get(topOfList - 1)) + ",");
            output.append(formatRGBColor(sortedColorPalette.get(topOfList - 2)) + "\n");

            writer.write(output.toString());
        }

        writer.close();
    }

    private static ArrayList<VBox> processImage(String urlString) throws IOException {
        BufferedImage image = downloadImage(urlString);

        CMap result = ImageAnalysis.getColorMap(image, 10);
        Collections.sort(result.vboxes, (x,y) -> x.compareTo(y));

        return result.vboxes;
    }

    private static BufferedImage downloadImage(String urlString) throws IOException {
        URL url = new URL(urlString);
        return ImageIO.read(url);
    }

    private static String formatRGBColor(VBox vbox) {
        int[] rgb = vbox.avg(false);

        // Create color String representations
        return createRGBHexString(rgb);
    }

    /**
     * Creates an HTML hex color code for the given RGB array (e.g. <code>#ff0000</code> for red).
     *
     * @param rgb
     *            the RGB array
     *
     * @return the HTML hex color code
     */
    private static String createRGBHexString(int[] rgb) {
        String rgbHex = Integer.toHexString(rgb[0] << 16 | rgb[1] << 8 | rgb[2]);

        // Left-pad with 0s
        int length = rgbHex.length();
        if (length < 6) {
            rgbHex = "00000".substring(0, 6 - length) + rgbHex;
        }

        return "#" + rgbHex.toUpperCase();
    }

}
