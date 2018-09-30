package mmcq.bulk.analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class MMCQ {

    private static final int SIGBITS = 5;                       // # of bits per color value (rgb). Losing 3 least sig
                                                                // bits for each RGB value has little effect on final color
    private static final int RSHIFT = 8 - SIGBITS;              // represents the shift amount to go from 8-bit input
                                                                // color value to what we expect (SIGBITS)
    private static final int MULT = 1 << RSHIFT;
    private static final int HISTOSIZE = 1 << (3 * SIGBITS);    // size of the histogram for RGB color space
                                                                // 3 = three colors (RGB)
    private static final int VBOX_LENGTH = 1 << SIGBITS;
    private static final double FRACT_BY_POPULATION = 0.75;     // fraction of VBoxes we divide based on
                                                                // population(# of pixels in VBox)
    private static final int MAX_ITERATIONS = 1000;

    /**
     * Get reduced-space color index for a pixel.
     * 
     * @param r
     *            the red value
     * @param g
     *            the green value
     * @param b
     *            the blue value
     * 
     * @return the color index
     */
    static int getColorIndex(int r, int g, int b) {
        return (r << (2 * SIGBITS)) + (g << SIGBITS) + b;
    }

    /**
     * 3D color space box. Each color value represents a side of the box.
     * Example: from point rMin -> point rMax
     */
    public static class VBox {
        int rMin;
        int rMax;
        int gMin;
        int gMax;
        int bMin;
        int bMax;

        private final int[] histo;  // histogram of color values and their counts from this VBox

        private int[] _avg;
        private Integer _volume;    // calculated volume of our rbg color space box
        private Integer _count;     // number of pixels

        public VBox(int rMin, int rMax, int gMin, int gMax, int bMin, int bMax, int[] histo) {
            this.rMin = rMin;
            this.rMax = rMax;
            this.gMin = gMin;
            this.gMax = gMax;
            this.bMin = bMin;
            this.bMax = bMax;

            this.histo = histo;
        }

        @Override
        public String toString() {
            return "rMin: " + rMin + " / rMax: " + rMax + " / gMin: " + gMin + " / gMax: " + gMax + " / bMin: " + bMin
                    + " / bMax: " + bMax;
        }

        // volume of this VBox calculated by our rgb side values
        public int volume(boolean force) {
            if (_volume == null || force) {
                _volume = ((rMax - rMin + 1) * (gMax - gMin + 1) * (bMax - bMin + 1));
            }

            return _volume;
        }

        // number of pixels in this VBox
        public int count(boolean force) {
            if (_count == null || force) {
                int numOfPixels = 0;
                int index;

                for (int r = rMin; r <= rMax; r++) {
                    for (int g = gMin; g <= gMax; g++) {
                        for (int b = bMin; b <= bMax; b++) {
                            index = getColorIndex(r, g, b);

                            // get the pixel count for this rgb value from the histogram
                            numOfPixels += histo[index];
                        }
                    }
                }

                _count = numOfPixels;
            }

            return _count;
        }

        @Override
        public VBox clone() {
            return new VBox(rMin, rMax, gMin, gMax, bMin, bMax, histo);
        }

        // finds the average color represented by the VBox
        public int[] avg(boolean force) {
            if (_avg == null || force) {
                int totalPixels = 0;

                int rsum = 0;
                int gsum = 0;
                int bsum = 0;

                int histoVal, histoIndex;

                for (int r = rMin; r <= rMax; r++) {
                    for (int g = gMin; g <= gMax; g++) {
                        for (int b = bMin; b <= bMax; b++) {
                            histoIndex = getColorIndex(r, g, b);
                            histoVal = histo[histoIndex];
                            totalPixels += histoVal;
                            rsum += (histoVal * (r + 0.5) * MULT);
                            gsum += (histoVal * (g + 0.5) * MULT);
                            bsum += (histoVal * (b + 0.5) * MULT);
                        }
                    }
                }

                if (totalPixels > 0) {
                    _avg = new int[] {rsum/totalPixels,
                            gsum/totalPixels,
                            bsum/totalPixels};
                } else {
                    _avg = new int[] {MULT * (rMin + rMax + 1) / 2,
                            MULT * (gMin + gMax + 1) / 2,
                            MULT * (bMin + bMax + 1) / 2};
                }
            }

            return _avg;
        }

        // Is the input color represented in the VBox?
        public boolean contains(int[] pixel) {
            int rval = pixel[0] >> RSHIFT;
            int gval = pixel[1] >> RSHIFT;
            int bval = pixel[2] >> RSHIFT;

            return (rval >= rMin && rval <= rMax && gval >= gMin && gval <= gMax && bval >= bMin
                    && bval <= bMax);
        }
    }

    /**
     * Color map.
     */
    public static class CMap {

        public final ArrayList<VBox> vboxes = new ArrayList<>();

        public void push(VBox box) {
            vboxes.add(box);
        }

        // gets the average color from each VBox
        public int[][] palette() {
            int numVBoxes = vboxes.size();
            int[][] palette = new int[numVBoxes][];
            for (int i = 0; i < numVBoxes; i++) {
                palette[i] = vboxes.get(i).avg(false);
            }
            return palette;
        }

        public int size() {
            return vboxes.size();
        }

        // finds the closest color represented by the palette
        public int[] map(int[] color) {
            int numVBoxes = vboxes.size();
            for (int i = 0; i < numVBoxes; i++) {
                VBox vbox = vboxes.get(i);
                if (vbox.contains(color)) {
                    return vbox.avg(false);
                }
            }
            return nearest(color);
        }

        // try to find the average color from one of the Vboxes
        // that is closest to our input color
        public int[] nearest(int[] inputColor) {
            double minResult = Double.MAX_VALUE;
            double result;
            int[] resultColor = null;

            int numVBoxes = vboxes.size();
            for (int i = 0; i < numVBoxes; i++) {
                int[] vbColor = vboxes.get(i).avg(false);

                // determine how large the gap between the RGB values is
                // result equal to zero means the values are equal
                // lower result means the values are closer
                // larger result means the values are further apart
                result = Math.sqrt(
                        Math.pow(inputColor[0] - vbColor[0], 2) + Math.pow(inputColor[1] - vbColor[1], 2)
                                + Math.pow(inputColor[2] - vbColor[2], 2));
                // trying to minimize the gap between RGB values
                if (result < minResult) {
                    minResult = result;
                    resultColor = vbColor;
                }
            }
            return resultColor;
        }

    }

    /**
     * Histo (1-d array, giving the number of pixels in each quantized region of color space), or
     * null on error.
     */
    private static int[] getHisto(int[][] pixels) {
        int[] histo = new int[HISTOSIZE];
        int index, rval, gval, bval;

        int numPixels = pixels.length;
        for (int i = 0; i < numPixels; i++) {
            int[] pixel = pixels[i];

            // we are receiving 8-bit color values, so we need to shift to meet
            // bit representation expectations set by SIGBITS
            rval = pixel[0] >> RSHIFT;
            gval = pixel[1] >> RSHIFT;
            bval = pixel[2] >> RSHIFT;
            index = getColorIndex(rval, gval, bval);
            histo[index]++;
        }
        return histo;
    }

    private static VBox vboxFromPixels(int[][] pixels, int[] histo) {
        int rmin = 1000000, rmax = 0;
        int gmin = 1000000, gmax = 0;
        int bmin = 1000000, bmax = 0;

        int rval, gval, bval;

        // find min/max
        int numPixels = pixels.length;
        for (int i = 0; i < numPixels; i++) {
            int[] pixel = pixels[i];
            rval = pixel[0] >> RSHIFT;
            gval = pixel[1] >> RSHIFT;
            bval = pixel[2] >> RSHIFT;

            if (rval < rmin) {
                rmin = rval;
            } else if (rval > rmax) {
                rmax = rval;
            }

            if (gval < gmin) {
                gmin = gval;
            } else if (gval > gmax) {
                gmax = gval;
            }

            if (bval < bmin) {
                bmin = bval;
            } else if (bval > bmax) {
                bmax = bval;
            }
        }

        return new VBox(rmin, rmax, gmin, gmax, bmin, bmax, histo);
    }

    private static VBox[] medianCutApply(int[] histo, VBox vbox) {
        if (vbox.count(false) == 0) {
            return null;
        }

        // only one pixel, no split
        if (vbox.count(false) == 1) {
            return new VBox[] {vbox.clone(), null};
        }

        int rw = vbox.rMax - vbox.rMin + 1;
        int gw = vbox.gMax - vbox.gMin + 1;
        int bw = vbox.bMax - vbox.bMin + 1;
        int maxw = Math.max(Math.max(rw, gw), bw);

        // Find the partial sum arrays along the selected axis.
        int total = 0;
        int[] partialsum = new int[VBOX_LENGTH];
        Arrays.fill(partialsum, -1); // -1 = not set / 0 = 0
        int[] lookaheadsum = new int[VBOX_LENGTH];
        Arrays.fill(lookaheadsum, -1); // -1 = not set / 0 = 0
        int i, j, k, sum, index;

        if (maxw == rw) {
            for (i = vbox.rMin; i <= vbox.rMax; i++) {
                sum = 0;
                for (j = vbox.gMin; j <= vbox.gMax; j++) {
                    for (k = vbox.bMin; k <= vbox.bMax; k++) {
                        index = getColorIndex(i, j, k);
                        sum += histo[index];
                    }
                }
                total += sum;
                partialsum[i] = total;
            }
        } else if (maxw == gw) {
            for (i = vbox.gMin; i <= vbox.gMax; i++) {
                sum = 0;
                for (j = vbox.rMin; j <= vbox.rMax; j++) {
                    for (k = vbox.bMin; k <= vbox.bMax; k++) {
                        index = getColorIndex(j, i, k);
                        sum += histo[index];
                    }
                }
                total += sum;
                partialsum[i] = total;
            }
        } else
        /* maxw == bw */
        {
            for (i = vbox.bMin; i <= vbox.bMax; i++) {
                sum = 0;
                for (j = vbox.rMin; j <= vbox.rMax; j++) {
                    for (k = vbox.gMin; k <= vbox.gMax; k++) {
                        index = getColorIndex(j, k, i);
                        sum += histo[index];
                    }
                }
                total += sum;
                partialsum[i] = total;
            }
        }

        for (i = 0; i < VBOX_LENGTH; i++) {
            if (partialsum[i] != -1) {
                lookaheadsum[i] = total - partialsum[i];
            }
        }

        // determine the cut planes
        return maxw == rw ? doCut('r', vbox, partialsum, lookaheadsum, total)
                : maxw == gw ? doCut('g', vbox, partialsum, lookaheadsum, total)
                        : doCut('b', vbox, partialsum, lookaheadsum, total);
    }

    private static VBox[] doCut(
            char color,
            VBox vbox,
            int[] partialsum,
            int[] lookaheadsum,
            int total) {
        int vbox_dim1;
        int vbox_dim2;

        if (color == 'r') {
            vbox_dim1 = vbox.rMin;
            vbox_dim2 = vbox.rMax;
        } else if (color == 'g') {
            vbox_dim1 = vbox.gMin;
            vbox_dim2 = vbox.gMax;
        } else
        /* color == 'b' */
        {
            vbox_dim1 = vbox.bMin;
            vbox_dim2 = vbox.bMax;
        }

        int left, right;
        VBox vbox1 = null, vbox2 = null;
        int d2, count2;

        for (int i = vbox_dim1; i <= vbox_dim2; i++) {
            if (partialsum[i] > total / 2) {
                vbox1 = vbox.clone();
                vbox2 = vbox.clone();

                left = i - vbox_dim1;
                right = vbox_dim2 - i;

                if (left <= right) {
                    d2 = Math.min(vbox_dim2 - 1, ~~(i + right / 2));
                } else {
                    // 2.0 and cast to int is necessary to have the same behaviour as in JavaScript
                    d2 = Math.max(vbox_dim1, ~~((int) (i - 1 - left / 2.0)));
                }

                // avoid 0-count boxes
                while (d2 < 0 || partialsum[d2] <= 0) {
                    d2++;
                }
                count2 = lookaheadsum[d2];
                while (count2 == 0 && d2 > 0 && partialsum[d2 - 1] > 0) {
                    count2 = lookaheadsum[--d2];
                }

                // set dimensions
                if (color == 'r') {
                    vbox1.rMax = d2;
                    vbox2.rMin = d2 + 1;
                } else if (color == 'g') {
                    vbox1.gMax = d2;
                    vbox2.gMin = d2 + 1;
                } else
                /* color == 'b' */
                {
                    vbox1.bMax = d2;
                    vbox2.bMin = d2 + 1;
                }

                return new VBox[] {vbox1, vbox2};
            }
        }

        throw new RuntimeException("VBox can't be cut");
    }

    // maxcolors: the color amount we want in the palette
    public static CMap quantize(int[][] pixels, int maxcolors) {
        // short-circuit
        if (pixels.length == 0 || maxcolors < 2 || maxcolors > 256) {
            return null;
        }

        int[] histo = getHisto(pixels);

        // get the beginning vbox from the colors
        VBox vbox = vboxFromPixels(pixels, histo);
        ArrayList<VBox> pq = new ArrayList<>();
        pq.add(vbox);

        // Round up to have the same behaviour as in JavaScript
        int target = (int) Math.ceil(FRACT_BY_POPULATION * maxcolors);

        // first set of colors, sorted by population
        iter(pq, COMPARATOR_COUNT, target, histo);

        // Re-sort by the product of pixel occupancy times the size in color space.
        Collections.sort(pq, COMPARATOR_PRODUCT);

        // next set - generate the median cuts using the (npix * vol) sorting.
        iter(pq, COMPARATOR_PRODUCT, maxcolors - pq.size(), histo);

        // Reverse to put the highest elements first into the color map
        Collections.reverse(pq);

        // calculate the actual colors
        CMap cmap = new CMap();
        for (VBox vb : pq) {
            cmap.push(vb);
        }

        return cmap;
    }

    /**
     * Inner function to do the iteration.
     */
    private static void iter(List<VBox> orgVboxes, Comparator<VBox> comparator, int target, int[] histo) {
        int ncolors = 1;
        int iterationCount = 0;
        VBox vbox;

        while (iterationCount < MAX_ITERATIONS) {
            vbox = orgVboxes.get(orgVboxes.size() - 1);
            if (vbox.count(false) == 0) {
                Collections.sort(orgVboxes, comparator);
                iterationCount++;
                continue;
            }
            orgVboxes.remove(orgVboxes.size() - 1);

            // do the cut
            VBox[] medianCutVboxes = medianCutApply(histo, vbox);
            VBox mcVbox1 = medianCutVboxes[0];
            VBox mcVbox2 = medianCutVboxes[1];

            if (mcVbox1 == null) {
                throw new RuntimeException("mcVbox1 not defined; shouldn't happen!");
            }

            orgVboxes.add(mcVbox1);
            if (mcVbox2 != null) {
                orgVboxes.add(mcVbox2);
                ncolors++;
            }
            Collections.sort(orgVboxes, comparator);

            if (ncolors >= target) {
                return;
            }
            if (iterationCount++ > MAX_ITERATIONS) {
                return;
            }
        }
    }

    public static final Comparator<VBox> COMPARATOR_COUNT = new Comparator<VBox>() {
        @Override
        public int compare(VBox a, VBox b) {
            return a.count(false) - b.count(false);
        }
    };

    private static final Comparator<VBox> COMPARATOR_PRODUCT = new Comparator<VBox>() {
        @Override
        public int compare(VBox a, VBox b) {
            int aCount = a.count(false);
            int bCount = b.count(false);
            int aVolume = a.volume(false);
            int bVolume = b.volume(false);

            // If count is 0 for both (or the same), sort by volume
            if (aCount == bCount) {
                return aVolume - bVolume;
            }

            // Otherwise sort by products
            return Long.compare((long) aCount * aVolume, (long) bCount * bVolume);
        }
    };

}
