package mmcq.bulk.analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class MMCQ {

    /**
     * what we've defined as the standard intensity value representation for most hardware
     */
    private static final int HARDWARE_INTENSITY_NUM_OF_BITS = 8;
    /**
     * # of bits per intensity value (rgb). Losing 3 least sigbits for each RGB value has little effect on final color
     */
    private static final int SIGBITS = 5;
    /**
     * represents the shift amount to go from the hardware intensity bit representation to what we expect (SIGBITS)
     */
    private static final int RSHIFT = HARDWARE_INTENSITY_NUM_OF_BITS - SIGBITS;
    /**
     * helps convert intensity to value of the range represented for HARDWARE_INTENSITY_NUM_OF_BITS to drive hardware pixel
     *
     * EXAMPLE: 8-bit hardware representation, our intensity value is 5-bits.
     * 5-bit value is represented at each increment of 5 from 0-255
     * Our step value is 2^3 so that intensity=5 is actually 5*(2^3) so it's correctly represented in the 0-255 hardware scale.
     */
    private static final int STEP_VALUE = 1 << RSHIFT;
    /**
     * size of the histogram for RGB color space 3 = three colors (RGB)
     */
    private static final int HISTOSIZE = 1 << (3 * SIGBITS);
    /**
     * number of possible intensity values for red or blue or green axis
     */
    private static final int VBOX_LENGTH = 1 << SIGBITS;
    /**
     * fraction of VBoxes we divide based on population(# of pixels in VBox)
     */
    private static final double FRACT_BY_POPULATION = 0.75;
    private static final int MAX_ITERATIONS = 1000;
    /**
     * when we do max-min subtraction, a color represented is lost so we need to add it back in
     */
    private static final int ADD_LOST_COLOR = 1;
    /**
     * when dividing by 2, force result to round up
     */
    private static final int FORCE_ROUND_UP = 1;

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
     * 3D color space box. Each color value represents the axises of the box.
     * Example: from point rMin -> point rMax is the red axis
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

        /**
         * volume of this VBox calculated by our rgb axis values
         */
        public int volume(boolean force) {
            if (_volume == null || force) {
                _volume = ((rMax - rMin + ADD_LOST_COLOR) * (gMax - gMin + ADD_LOST_COLOR) * (bMax - bMin + ADD_LOST_COLOR));
            }

            return _volume;
        }

        /**
         * Number of pixels in this VBox
         * @param force
         *              force a re-count
         */
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

        /**
         * Finds the average color represented by the VBox
         * @param force
         *              force getting the average even if it exists
         */
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
                            // convert value intensity to value in range 0-255 to drive hardware pixel (8-bit)
                            // Example: 5-bit value is represented at each increment of 5 from 0-255
                            // 0.5 is to round to next intensity to avoid truncation when we divide against
                            // total pixel amount later (which would lessen the intensity of the value)
                            rsum += (histoVal * ((r + 0.5) * STEP_VALUE));
                            gsum += (histoVal * ((g + 0.5) * STEP_VALUE));
                            bsum += (histoVal * ((b + 0.5) * STEP_VALUE));
                        }
                    }
                }

                if (totalPixels > 0) {
                    _avg = new int[] {rsum/totalPixels,
                            gsum/totalPixels,
                            bsum/totalPixels};
                } else {
                    _avg = new int[] {STEP_VALUE * (rMin + rMax + FORCE_ROUND_UP) / 2,
                            STEP_VALUE * (gMin + gMax + FORCE_ROUND_UP) / 2,
                            STEP_VALUE * (bMin + bMax + FORCE_ROUND_UP) / 2};
                }
            }

            return _avg;
        }

        /**
         * Is the input color represented in the VBox?
         * @param pixel
         */
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

        /**
         * Gets the average color from each VBox
         */
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

        /**
         * Finds the closest color represented by the palette
         * @param color
         */
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

        /**
         * Try to find the average color from one of the Vboxes that is closest to our input color
         */
        public int[] nearest(int[] inputColor) {
            double minResult = Double.MAX_VALUE;
            double result;
            int[] resultColor = null;

            int numVBoxes = vboxes.size();
            for (int i = 0; i < numVBoxes; i++) {
                int[] vbColor = vboxes.get(i).avg(false);

                /**
                 * determine how large the gap between the RGB values is
                 * result = 0; the values are equal
                 * lower result means the values are closer
                 * higher result means the values are further apart
                 */
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

            /**
             * we are receiving 8-bit color values, so we need to shift to meet
             * bit representation expectations set by SIGBITS
             */
            rval = pixel[0] >> RSHIFT;
            gval = pixel[1] >> RSHIFT;
            bval = pixel[2] >> RSHIFT;
            index = getColorIndex(rval, gval, bval);
            histo[index]++;
        }
        return histo;
    }

    private static VBox vboxFromPixels(int[][] pixels, int[] histo) {
        int rMin = 1000000, rMax = 0;
        int gMin = 1000000, gMax = 0;
        int bMin = 1000000, bMax = 0;

        int rVal, gVal, bVal;

        /**
         * find min/max to determine the range of the VBox we are creating
         * aka the points to measure the axises of our color space (ex. red axis: rMin->rMax)
         */
        int numPixels = pixels.length;
        for (int i = 0; i < numPixels; i++) {
            int[] pixel = pixels[i];
            /**
             * we are receiving 8-bit color values, so we need to shift to meet
             * bit representation expectations set by SIGBITS
             */
            rVal = pixel[0] >> RSHIFT;
            gVal = pixel[1] >> RSHIFT;
            bVal = pixel[2] >> RSHIFT;

            if (rVal < rMin) {
                rMin = rVal;
            } else if (rVal > rMax) {
                rMax = rVal;
            }

            if (gVal < gMin) {
                gMin = gVal;
            } else if (gVal > gMax) {
                gMax = gVal;
            }

            if (bVal < bMin) {
                bMin = bVal;
            } else if (bVal > bMax) {
                bMax = bVal;
            }
        }

        return new VBox(rMin, rMax, gMin, gMax, bMin, bMax, histo);
    }

    private static VBox[] medianCutApply(int[] histo, VBox orgVbox) {
        if (orgVbox.count(false) == 0) {
            return null;
        }

        // only one pixel, no split. return original vbox
        if (orgVbox.count(false) == 1) {
            return new VBox[] {orgVbox.clone(), null};
        }

        // calculate axis values of the VBox (color space)
        int rAxis = orgVbox.rMax - orgVbox.rMin + ADD_LOST_COLOR;
        int gAxis = orgVbox.gMax - orgVbox.gMin + ADD_LOST_COLOR;
        int bAxis = orgVbox.bMax - orgVbox.bMin + ADD_LOST_COLOR;
        int maxAxis = Math.max(Math.max(rAxis, gAxis), bAxis);

        // Find the partial sum arrays along the selected axis.
        int total = 0;

        /**
         * the sum of the previous(if exists) and the pixel count of each present index of range(longest axis)
         */
        int[] partialsum = new int[VBOX_LENGTH];
        Arrays.fill(partialsum, -1); // -1 = not set / 0 = 0

        int[] lookaheadsum = new int[VBOX_LENGTH];
        Arrays.fill(lookaheadsum, -1); // -1 = not set / 0 = 0
        int sum, index;

        // calculate cumulative count of pixels for groups of hues ordered by range of longest axis
        if (maxAxis == rAxis) {
            for (int r = orgVbox.rMin; r <= orgVbox.rMax; r++) {
                sum = 0;
                for (int g = orgVbox.gMin; g <= orgVbox.gMax; g++) {
                    for (int b = orgVbox.bMin; b <= orgVbox.bMax; b++) {
                        index = getColorIndex(r, g, b);
                        sum += histo[index];
                    }
                }
                total += sum;
                partialsum[r] = total;
            }
        } else if (maxAxis == gAxis) {
            for (int g = orgVbox.gMin; g <= orgVbox.gMax; g++) {
                sum = 0;
                for (int r = orgVbox.rMin; r <= orgVbox.rMax; r++) {
                    for (int b = orgVbox.bMin; b <= orgVbox.bMax; b++) {
                        index = getColorIndex(r, g, b);
                        sum += histo[index];
                    }
                }
                total += sum;
                partialsum[g] = total;
            }
        } else
        /* maxAxis == bAxis */
        {
            for (int b = orgVbox.bMin; b <= orgVbox.bMax; b++) {
                sum = 0;
                for (int r = orgVbox.rMin; r <= orgVbox.rMax; r++) {
                    for (int g = orgVbox.gMin; g <= orgVbox.gMax; g++) {
                        index = getColorIndex(r, g, b);
                        sum += histo[index];
                    }
                }
                total += sum;
                partialsum[b] = total;
            }
        }

        for (int i = 0; i < VBOX_LENGTH; i++) {
            if (partialsum[i] != -1) {
                /**
                 * the pixel sum of the groups of hues in future indices of range(longest axis)
                 */
                lookaheadsum[i] = total - partialsum[i];
            }
        }

        // determine the cut planes
        return maxAxis == rAxis ? doCut('r', orgVbox, partialsum, lookaheadsum, total)
                : maxAxis == gAxis ? doCut('g', orgVbox, partialsum, lookaheadsum, total)
                        : doCut('b', orgVbox, partialsum, lookaheadsum, total);
    }

    private static VBox[] doCut(
            char color,
            VBox orgVbox,
            int[] partialsum,
            int[] lookaheadsum,
            int total) {
        int vbox_min;
        int vbox_max;

        if (color == 'r') {
            vbox_min = orgVbox.rMin;
            vbox_max = orgVbox.rMax;
        } else if (color == 'g') {
            vbox_min = orgVbox.gMin;
            vbox_max = orgVbox.gMax;
        } else
        /* color == 'b' */
        {
            vbox_min = orgVbox.bMin;
            vbox_max = orgVbox.bMax;
        }

        int left, right;
        VBox vbox1 = null, vbox2 = null;
        int vboxCutVal, count2;

        /**
         * find first value in range of min to max whose cumulative pixel count is greater than half
         * the total pixel count.
         * This insures that when we split the original vbox, the left(vbox1) has a greater population than right(vbox2)
         */
        for (int i = vbox_min; i <= vbox_max; i++) {
            if (partialsum[i] > total / 2) {
                vbox1 = orgVbox.clone();
                vbox2 = orgVbox.clone();

                left = i - vbox_min;
                right = vbox_max - i;

                if (left <= right) {
                    vboxCutVal = Math.min(vbox_max - 1, (i + (right / 2)));
                } else {
                    // 2.0 and cast to int is necessary to have the same behaviour as in JavaScript
                    vboxCutVal = Math.max(vbox_min, ((int) (i - 1 - (left / 2.0))));
                }

                // avoid 0-count boxes
                while (vboxCutVal < 0 || partialsum[vboxCutVal] <= 0) {
                    vboxCutVal++;
                }
                count2 = lookaheadsum[vboxCutVal];
                // count2 == 0: means count of pixels for value in range of current to max are empty
                // vboxCutVal > 0: means there exists a pixel count partial sum below our split point (vboxCutVal)
                // If needed, move split down to prevent vbox2 from containing no pixels
                while (count2 == 0 && vboxCutVal > 0 && partialsum[vboxCutVal - 1] > 0) {
                    count2 = lookaheadsum[--vboxCutVal];
                }

                // set dimensions
                if (color == 'r') {
                    vbox1.rMax = vboxCutVal;
                    vbox2.rMin = vboxCutVal + 1;
                } else if (color == 'g') {
                    vbox1.gMax = vboxCutVal;
                    vbox2.gMin = vboxCutVal + 1;
                } else
                /* color == 'b' */
                {
                    vbox1.bMax = vboxCutVal;
                    vbox2.bMin = vboxCutVal + 1;
                }

                return new VBox[] {vbox1, vbox2};
            }
        }

        throw new RuntimeException("VBox can't be cut");
    }

    /**
     *
     * @param pixels
     * @param maxcolors
     *                  the amount of colors we want in the result palette
     * @return
     */
    public static CMap quantize(int[][] pixels, int maxcolors) {
        // short-circuit
        if (pixels.length == 0 || maxcolors < 2 || maxcolors > 256) {
            return null;
        }

        int[] histo = getHisto(pixels);

        // get the beginning vbox from the colors
        VBox vbox = vboxFromPixels(pixels, histo);
        ArrayList<VBox> vboxList = new ArrayList<>();
        vboxList.add(vbox);

        // for a fraction of the color palette generation, we subdivide based on population in VBox
        // Round up to have the same behaviour as in JavaScript
        int fractOfTarget = (int) Math.ceil(FRACT_BY_POPULATION * maxcolors);

        // first set of colors, sorted by population
        iter(vboxList, COMPARATOR_COUNT, fractOfTarget, histo);

        // Re-sort by the product of pixel occupancy times the size in color space.
        Collections.sort(vboxList, COMPARATOR_PRODUCT);

        // next set (rest of color palette generation) - generate the median cuts using the (npix * vol) sorting.
        iter(vboxList, COMPARATOR_PRODUCT, maxcolors - vboxList.size(), histo);

        // Reverse to put the highest elements first into the color map
        Collections.reverse(vboxList);

        // calculate the actual colors
        CMap cmap = new CMap();
        for (VBox vb : vboxList) {
            cmap.push(vb);
        }

        return cmap;
    }

    /**
     * Inner function to do the iteration.
     */
    private static void iter(List<VBox> orgVboxes, Comparator<VBox> comparator, int target, int[] histo) {
        int numOfPaletteColors = 1;
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
                numOfPaletteColors++;
            }
            Collections.sort(orgVboxes, comparator);

            if (numOfPaletteColors >= target) {
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
