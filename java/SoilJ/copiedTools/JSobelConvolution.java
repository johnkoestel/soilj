package SoilJ.copiedTools;

/**
 * SoilJ.copiedTools is a collection of classes needed to run SoilJ, 
 * that were minimally modified 2014 - 2016 by John Koestel to to make it accessible from 
 * within SoilJ. 
 *
 * The original file is
 *  
 * <p>Title: Sobel Folding 3D Plugin for ImageJ</p>
 *
 * <p>Description: Convolves an image with a Sobel Kernel</p>
 *
 * <p>Copyright: Copyright (c) 2008</p>
 *
 * <p>Company: MPI-CBG</p>
 *
 * <p>License: GPL
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 2
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * @author Stephan Preibisch
 * @author John Koestel (refactoring)
 * @version 1.0
 */


import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

public class JSobelConvolution implements PlugIn
{
    public JSobelConvolution()
    {
    }

    private FloatArray img;

    /**
     * This method will be called when running the PlugIn, it coordinates the main process.
     *
     * @param args UNUSED

     * @author   Stephan Preibisch
     */
    public void run(String arg) {}
    /*{
        ImagePlus imp = WindowManager.getCurrentImage();

        if (null == imp)
        {
            IJ.log("No images open.");
            return;
        }

        //
        // Show the dialog
        //
        String methodList[] = {"Mirrored", "Black", "Minium Value", "Maximum Value", "Average Value", "Value Given Below"};

        GenericDialog gd = new GenericDialog("Sobel Convolution 3D");

        gd.addChoice("Behaviour at image borders, assume pixels outside the image to be ", methodList, methodList[0]);
        gd.addNumericField("Background value (only needed if this option is selected above)", 0, 3);

        gd.showDialog();

        if (gd.wasCanceled())
            return;

        //
        // check whether it is a 3D image and not RGB
        //
        ImageStack stack = imp.getStack();
        img = StackToFloatArray(stack);
        if (img == null)
        {
            IJ.error("RGB not supported yet. You can split into seperate channels and merge afterwards again...sorry!");
            return;
        }

        String behaviour = gd.getNextChoice();

        float value = 0;

        if (behaviour.equals(methodList[methodList.length - 1]))
            value = (float) gd.getNextNumber();
        else if (behaviour.equals("Minium Value"))
            value = getMin(img);
        else if (behaviour.equals("Maximum Value"))
            value = getMax(img);
        else if (behaviour.equals("Average Value"))
            value = getAverage(img);

        //
        // compute Gaussian
        //
        if (behaviour.equals(methodList[0]))
            img = computeSobelMirror((FloatArray3D) img);
        else
            img = computeSobelValue((FloatArray3D) img, value);

        FloatArrayToStack((FloatArray3D) img, "Sobel Filtered Image", getMin(img), getMax(img)).show();
    }*/
    
    public ImagePlus getMeWhatIWant(ImagePlus nowTiff) {
    	
    	ImagePlus outTiff = new ImagePlus();
    	
    	//apply the Sobel 3D filter with Black background and return image 
    	ImageStack stack = nowTiff.getStack();
        img = StackToFloatArray(stack);
        
        img = computeSobelMirror((FloatArray3D) img);
    	
        outTiff = FloatArrayToStack((FloatArray3D) img, "Sobel Filtered Image", getMin(img), getMax(img));
                
        return outTiff;
    }

    private static float getAverage(FloatArray img)
    {
        double avg = 0;

        for (int i = 0; i < img.data.length; i++)
            avg += img.data[i];

        return (float) (avg / (double) img.data.length);
    }

    private static float getMin(FloatArray img)
    {
        float min = Float.MAX_VALUE;

        for (int i = 0; i < img.data.length; i++)
        {
            float value = img.data[i];
            if (value < min)
                min = value;
        }

        return min;
    }

    private static float getMax(FloatArray img)
    {
        float max = Float.MIN_VALUE;

        for (int i = 0; i < img.data.length; i++)
        {
            float value = img.data[i];
            if (value > max)
                max = value;
        }

        return max;
    }

    /**
     * This method converts my FloatArray2D to an ImageJ ImagePlus
     *
     * @param image The image as FloatArray2D
     * @param name The name of the ImagePlus
     * @param min Lowest brightness value that will be displayed (see Brightness&Contrast in Imagej)
     * @param max Highest brightness value that will be displayed (set both to zero for automatic)
     * @return ImagePlus The ImageJ image
     *
     * @author   Stephan Preibisch
     */
    public static ImagePlus FloatArrayToImagePlus(FloatArray2D image, String name, float min, float max)
    {
        ImagePlus imp = IJ.createImage(name, "32-Bit Black", image.width, image.height, 1);
        FloatProcessor ip = (FloatProcessor) imp.getProcessor();
        FloatArrayToFloatProcessor(ip, image);

        if (min == max)
            ip.resetMinAndMax();
        else
            ip.setMinAndMax(min, max);

        imp.updateAndDraw();

        return imp;
    }

    /**
     * This method converts my FloatArray2D to an ImageJ ImageProcessor
     *
     * @param ImageProcessor Will be overwritten with the img from the FloatArray2D
     * @param FloatArray2D The image as FloatArray2D
     * @return
     *
     * @author   Stephan Preibisch
     */
    public static void FloatArrayToFloatProcessor(ImageProcessor ip, FloatArray2D pixels)
    {
        float[] img = new float[pixels.width * pixels.height];

        int count = 0;
        for (int y = 0; y < pixels.height; y++)
            for (int x = 0; x < pixels.width; x++)
                img[count] = pixels.data[count++];

        ip.setPixels(img);
        ip.resetMinAndMax();
    }


    /**
     * This method converts my FloatArray3D to an ImageJ image stack packed into an ImagePlus
     *
     * @param image The image as FloatArray3D
     * @param name The name of the ImagePlus
     * @param min Lowest brightness value that will be displayed (see Brightness&Contrast in Imagej)
     * @param max Highest brightness value that will be displayed (set both to zero for automatic)
     * @return ImagePlus The ImageJ image
     *
     * @author   Stephan Preibisch
     */
    public ImagePlus FloatArrayToStack(FloatArray3D image, String name, float min, float max)
    {
        int width = image.width;
        int height = image.height;
        int nstacks = image.depth;

        ImageStack stack = new ImageStack(width, height);

        for (int slice = 0; slice < nstacks; slice++)
        {
            ImagePlus impResult = IJ.createImage(name, "32-Bit Black", width, height, 1);
            ImageProcessor ipResult = impResult.getProcessor();
            float[] sliceImg = new float[width * height];

            for (int x = 0; x < width; x++)
                for (int y = 0; y < height; y++)
                    sliceImg[y * width + x] = image.get(x, y, slice);

            ipResult.setPixels(sliceImg);

            if (min == max)
                ipResult.resetMinAndMax();
            else
                ipResult.setMinAndMax(min, max);

            stack.addSlice("Slice " + slice, ipResult);
        }

        return new ImagePlus(name, stack);
    }
    
    private FloatArray3D computeSobelMirror(FloatArray3D input)
    {
            FloatArray3D output = new FloatArray3D(input.width, input.height, input.depth);

            float sobelX, sobelY, sobelZ;
            float v1, v2, v3, v4, v5, v6, v7, v8, v9, v10;

            for (int z = 0; z < input.depth; z++)
                    for (int y = 0; y < input.height; y++)
                            for (int x = 0; x < input.width; x++)
                            {
                                    // Sobel in X
                                    v1 = input.getMirror(x - 1, y, z - 1);
                                    v2 = input.getMirror(x + 1, y, z - 1);

                                    v3 = input.getMirror(x - 1, y - 1, z);
                                    v4 = input.getMirror(x - 1, y, z);
                                    v5 = input.getMirror(x - 1, y + 1, z);

                                    v6 = input.getMirror(x + 1, y - 1, z);
                                    v7 = input.getMirror(x + 1, y, z);
                                    v8 = input.getMirror(x + 1, y + 1, z);

                                    v9 = input.getMirror(x - 1, y, z + 1);
                                    v10 = input.getMirror(x + 1, y, z + 1);

                                    sobelX = v1 + ( -v2) +
                                                     v3 + (2 * v4) + v5 +
                                                     ( -v6) + ( -2 * v7) + ( -v8) +
                                                     v9 + ( -v10);

                                    // Sobel in Y
                                    v1 = input.getMirror(x, y - 1, z - 1);
                                    v2 = input.getMirror(x, y + 1, z - 1);

                                    v3 = input.getMirror(x - 1, y - 1, z);
                                    v4 = input.getMirror(x, y - 1, z);
                                    v5 = input.getMirror(x + 1, y - 1, z);

                                    v6 = input.getMirror(x - 1, y + 1, z);
                                    v7 = input.getMirror(x, y + 1, z);
                                    v8 = input.getMirror(x + 1, y + 1, z);

                                    v9 = input.getMirror(x, y - 1, z + 1);
                                    v10 = input.getMirror(x, y + 1, z + 1);

                                    sobelY = v1 + ( -v2) +
                                                     v3 + (2 * v4) + v5 +
                                                     ( -v6) + ( -2 * v7) + ( -v8) +
                                                     v9 + ( -v10);

                                    // Sobel in Z
                                    v1 = input.getMirror(x, y - 1, z - 1);
                                    v2 = input.getMirror(x - 1, y, z - 1);
                                    v3 = input.getMirror(x, y, z - 1);
                                    v4 = input.getMirror(x + 1, y, z - 1);
                                    v5 = input.getMirror(x, y + 1, z - 1);

                                    v6 = input.getMirror(x, y - 1, z + 1);
                                    v7 = input.getMirror(x - 1, y, z + 1);
                                    v8 = input.getMirror(x, y, z + 1);
                                    v9 = input.getMirror(x + 1, y, z + 1);
                                    v10 = input.getMirror(x, y + 1, z + 1);

                                    sobelZ = v1 +
                                                     v2 + (2 * v3) + v4 +
                                                     v5 +
                                                     ( -v6) +
                                                     ( -v7) + ( -2 * v8) + ( -v9) +
                                                     ( -v10);

                                    output.set((float) Math.sqrt(Math.pow(sobelX, 2) + Math.pow(sobelY, 2) + Math.pow(sobelZ, 2)), x, y, z);
                            }

            return output;
    }

    private FloatArray3D computeSobelValue(FloatArray3D input, float value)
    {
            FloatArray3D output = new FloatArray3D(input.width, input.height, input.depth);

            float sobelX, sobelY, sobelZ;
            float v1, v2, v3, v4, v5, v6, v7, v8, v9, v10;

            for (int z = 0; z < input.depth; z++)
                    for (int y = 0; y < input.height; y++)
                            for (int x = 0; x < input.width; x++)
                            {
                                    // Sobel in X
                                    v1 = input.getValueOutOfImage(x - 1, y, z - 1, value);
                                    v2 = input.getValueOutOfImage(x + 1, y, z - 1, value);

                                    v3 = input.getValueOutOfImage(x - 1, y - 1, z, value);
                                    v4 = input.getValueOutOfImage(x - 1, y, z, value);
                                    v5 = input.getValueOutOfImage(x - 1, y + 1, z, value);

                                    v6 = input.getValueOutOfImage(x + 1, y - 1, z, value);
                                    v7 = input.getValueOutOfImage(x + 1, y, z, value);
                                    v8 = input.getValueOutOfImage(x + 1, y + 1, z, value);

                                    v9 = input.getValueOutOfImage(x - 1, y, z + 1, value);
                                    v10 = input.getValueOutOfImage(x + 1, y, z + 1, value);

                                    sobelX = v1 + ( -v2) +
                                                     v3 + (2 * v4) + v5 +
                                                     ( -v6) + ( -2 * v7) + ( -v8) +
                                                     v9 + ( -v10);

                                    // Sobel in Y
                                    v1 = input.getValueOutOfImage(x, y - 1, z - 1, value);
                                    v2 = input.getValueOutOfImage(x, y + 1, z - 1, value);

                                    v3 = input.getValueOutOfImage(x - 1, y - 1, z, value);
                                    v4 = input.getValueOutOfImage(x, y - 1, z, value);
                                    v5 = input.getValueOutOfImage(x + 1, y - 1, z, value);

                                    v6 = input.getValueOutOfImage(x - 1, y + 1, z, value);
                                    v7 = input.getValueOutOfImage(x, y + 1, z, value);
                                    v8 = input.getValueOutOfImage(x + 1, y + 1, z, value);

                                    v9 = input.getValueOutOfImage(x, y - 1, z + 1, value);
                                    v10 = input.getValueOutOfImage(x, y + 1, z + 1, value);

                                    sobelY = v1 + ( -v2) +
                                                     v3 + (2 * v4) + v5 +
                                                     ( -v6) + ( -2 * v7) + ( -v8) +
                                                     v9 + ( -v10);

                                    // Sobel in Z
                                    v1 = input.getValueOutOfImage(x, y - 1, z - 1, value);
                                    v2 = input.getValueOutOfImage(x - 1, y, z - 1, value);
                                    v3 = input.getValueOutOfImage(x, y, z - 1, value);
                                    v4 = input.getValueOutOfImage(x + 1, y, z - 1, value);
                                    v5 = input.getValueOutOfImage(x, y + 1, z - 1, value);

                                    v6 = input.getValueOutOfImage(x, y - 1, z + 1, value);
                                    v7 = input.getValueOutOfImage(x - 1, y, z + 1, value);
                                    v8 = input.getValueOutOfImage(x, y, z + 1, value);
                                    v9 = input.getValueOutOfImage(x + 1, y, z + 1, value);
                                    v10 = input.getValueOutOfImage(x, y + 1, z + 1, value);

                                    sobelZ = v1 +
                                                     v2 + (2 * v3) + v4 +
                                                     v5 +
                                                     ( -v6) +
                                                     ( -v7) + ( -2 * v8) + ( -v9) +
                                                     ( -v10);

                                    output.set((float) Math.sqrt(Math.pow(sobelX, 2) + Math.pow(sobelY, 2) + Math.pow(sobelZ, 2)), x, y, z);
                            }

            return output;
    }

    /**
     * This method convertes an ImageJ image stack to my FloatArray3D,
     * which is a one dimensional structure with methods for 3D access
     *
     * @param stack ImageJ image stack
     * @return FloatArray3D The image packed into a FloatArray3D
     *
     * @author   Stephan Preibisch
     */
    public FloatArray3D StackToFloatArray(ImageStack stack)
    {
        Object[] imageStack = stack.getImageArray();
        int width = stack.getWidth();
        int height = stack.getHeight();
        int nstacks = stack.getSize();

        if (imageStack == null || imageStack.length == 0)
        {
            IJ.log("Image Stack is empty.");
            return null;
        }

        if (imageStack[0] instanceof int[])
        {
            IJ.log("RGB images supported at the moment.");
            return null;
        }

        FloatArray3D pixels = new FloatArray3D(width, height, nstacks);
        //float[][][] pixels = new float[width][height][nstacks];
        int count;

        if (imageStack[0] instanceof byte[])
            for (int countSlice = 0; countSlice < nstacks; countSlice++)
            {
                byte[] pixelTmp = (byte[]) imageStack[countSlice];
                count = 0;

                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        pixels.data[pixels.getPos(x, y, countSlice)] = (float) (pixelTmp[count++] & 0xff);
            }
        else if (imageStack[0] instanceof short[])
            for (int countSlice = 0; countSlice < nstacks; countSlice++)
            {
                short[] pixelTmp = (short[]) imageStack[countSlice];
                count = 0;

                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        pixels.data[pixels.getPos(x, y, countSlice)] = (float) (pixelTmp[count++] & 0xffff);
            }
        else // instance of float[]
            for (int countSlice = 0; countSlice < nstacks; countSlice++)
            {
                float[] pixelTmp = (float[]) imageStack[countSlice];
                count = 0;

                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        pixels.data[pixels.getPos(x, y, countSlice)] = pixelTmp[count++];
            }

        return pixels;
    }

    /**
     * This method convertes an ImageJ ImageProcessor to my FloatArray2D,
     * which is a one dimensional structure with methods for 2D access
     *
     * @param stack ImageJ ImageProcessor
     * @return FloatArray2D The image packed into a FloatArray2D
     *
     * @author   Stephan Preibisch
     */
    public FloatArray2D ImageToFloatArray(ImageProcessor ip)
    {
        FloatArray2D image;
        Object pixelArray = ip.getPixels();
        int count = 0;

        if (ip instanceof ByteProcessor)
        {
            image = new FloatArray2D(ip.getWidth(), ip.getHeight());
            byte[] pixels = (byte[]) pixelArray;

            for (int y = 0; y < ip.getHeight(); y++)
                for (int x = 0; x < ip.getWidth(); x++)
                    image.data[count] = pixels[count++] & 0xff;
        }
        else if (ip instanceof ShortProcessor)
        {
            image = new FloatArray2D(ip.getWidth(), ip.getHeight());
            short[] pixels = (short[]) pixelArray;

            for (int y = 0; y < ip.getHeight(); y++)
                for (int x = 0; x < ip.getWidth(); x++)
                    image.data[count] = pixels[count++] & 0xffff;
        }
        else if (ip instanceof FloatProcessor)
        {
            image = new FloatArray2D(ip.getWidth(), ip.getHeight());
            float[] pixels = (float[]) pixelArray;

            for (int y = 0; y < ip.getHeight(); y++)
                for (int x = 0; x < ip.getWidth(); x++)
                    image.data[count] = pixels[count++];
        }
        else //RGB
        {
            IJ.log("RGB images not supported");
            image = null;
        }

        return image;
    }

    /**
     * This class is the abstract class for my FloatArrayXDs,
     * which are a one dimensional structures with methods for access in n dimensions
     *
     * @author   Stephan Preibisch
     */
    public abstract class FloatArray
    {
        public float data[] = null;
        public abstract FloatArray clone();
    }


    /**
     * The 2D implementation of the FloatArray
     *
     * @author   Stephan Preibisch
     */
    public class FloatArray2D extends FloatArray
    {
        public int width = 0;
        public int height = 0;

        private int doubleWidth, doubleHeight;

        public FloatArray2D(int width, int height)
        {
            data = new float[width * height];
            this.width = width;
            this.height = height;

            doubleWidth = 2 * width;
            doubleHeight = 2 * height;
        }

        public FloatArray2D(float[] data, int width, int height)
        {
            this.data = data;
            this.width = width;
            this.height = height;

            doubleWidth = 2 * width;
            doubleHeight = 2 * height;
        }

        public FloatArray2D clone()
        {
            FloatArray2D clone = new FloatArray2D(width, height);
            System.arraycopy(this.data, 0, clone.data, 0, this.data.length);
            return clone;
        }

        public int getPos(int x, int y)
        {
            return x + width * y;
        }

        public float get(int x, int y)
        {
            return data[getPos(x, y)];
        }

        public float getValueOutOfImage(int x, int y, float value)
        {
            if (x > 0 && x < width && y > 0 && y < height)
                return data[getPos(x, y)];
            else
                return value;
        }

        public float getFlipInRange(int x, int y)
        {
            if (x < 0) x = doubleWidth + x % doubleWidth;
            if (x >= doubleWidth) x = x % doubleWidth;
            if (x >= width) x = width - x % width - 1;

            if (y < 0) y = doubleHeight + y % doubleHeight;
            if (y >= doubleHeight) y = y % doubleHeight;
            if (y >= height) y = height - y % height - 1;

            return data[getPos(x, y)];
            //return get(flipInRange(x, width), flipInRange(y, height) );
        }

        public float getMirror(int x, int y)
        {
            if (x >= width)
                x = width - (x - width + 2);

            if (y >= height)
                y = height - (y - height + 2);

            if (x < 0)
            {
                int tmp = 0;
                int dir = 1;

                while (x < 0)
                {
                    tmp += dir;
                    if (tmp == width - 1 || tmp == 0)
                        dir *= -1;
                    x++;
                }
                x = tmp;
            }

            if (y < 0)
            {
                int tmp = 0;
                int dir = 1;

                while (y < 0)
                {
                    tmp += dir;
                    if (tmp == height - 1 || tmp == 0)
                        dir *= -1;
                    y++;
                }
                y = tmp;
            }

            return data[getPos(x, y)];
        }

        public float getZero(int x, int y)
        {
            if (x >= width)
                return 0;

            if (y >= height)
                return 0;

            if (x < 0)
                return 0;

            if (y < 0)
                return 0;

            return data[getPos(x, y)];
        }

        public void set(float value, int x, int y)
        {
            data[getPos(x, y)] = value;
        }
    }


    /**
     * The 3D implementation of the FloatArray
     *
     * @author   Stephan Preibisch
     */
    public class FloatArray3D extends FloatArray
    {
        public int width = 0;
        public int height = 0;
        public int depth = 0;

        private int doubleWidth, doubleHeight, doubleDepth;

        public FloatArray3D(float[] data, int width, int height, int depth)
        {
            this.data = data;
            this.width = width;
            this.height = height;
            this.depth = depth;

            doubleWidth = 2 * width;
            doubleHeight = 2 * height;
            doubleDepth = 2 * depth;
        }

        public FloatArray3D(int width, int height, int depth)
        {
            data = new float[width * height * depth];
            this.width = width;
            this.height = height;
            this.depth = depth;

            doubleWidth = 2 * width;
            doubleHeight = 2 * height;
            doubleDepth = 2 * depth;
        }

        public FloatArray3D clone()
        {
            FloatArray3D clone = new FloatArray3D(width, height, depth);
            System.arraycopy(this.data, 0, clone.data, 0, this.data.length);
            return clone;
        }

        public int getPos(int x, int y, int z)
        {
            return x + width * (y + z * height);
        }

        public float get(int x, int y, int z)
        {
            return data[getPos(x, y, z)];
        }

        public float getValueOutOfImage(int x, int y, int z, float value)
        {
            if (x > 0 && x < width && y > 0 && y < height && z > 0 && z < depth)
                return data[getPos(x, y, z)];
            else
                return value;
        }

        /**
         * Return the value at an arbitrary position, where the image data is flipped like that:
         *
         * Size = 3
         *
         * -4 -> 2
         * -3 -> 2
         * -2 -> 1
         * -1 -> 0
         * 0 -> 0
         * 1 -> 1
         * 2 -> 2
         * 3 -> 2
         * 4 -> 1
         * 5 -> 0
         * 6 -> 2
         *
         * @param x int x position
         * @param y int y position
         * @param z int z position
         * @return float The value
         */
        public float getFlipInRange(int x, int y, int z)
        {
            if (x < 0) x = doubleWidth + x % doubleWidth;
            if (x >= doubleWidth) x = x % doubleWidth;
            if (x >= width) x = width - x % width - 1;

            if (y < 0) y = doubleHeight + y % doubleHeight;
            if (y >= doubleHeight) y = y % doubleHeight;
            if (y >= height) y = height - y % height - 1;

            if (z < 0) z = doubleDepth + z % doubleDepth;
            if (z >= doubleDepth) z = z % doubleDepth;
            if (z >= depth) z = depth - z % depth - 1;

            return data[getPos(x, y, z)];
        }

        public float getMirror(int x, int y, int z)
        {
            if (x >= width)
                x = width - (x - width + 2);

            if (y >= height)
                y = height - (y - height + 2);

            if (z >= depth)
                z = depth - (z - depth + 2);

            if (x < 0)
            {
                int tmp = 0;
                int dir = 1;

                while (x < 0)
                {
                    tmp += dir;
                    if (tmp == width - 1 || tmp == 0)
                        dir *= -1;
                    x++;
                }
                x = tmp;
            }

            if (y < 0)
            {
                int tmp = 0;
                int dir = 1;

                while (y < 0)
                {
                    tmp += dir;
                    if (tmp == height - 1 || tmp == 0)
                        dir *= -1;
                    y++;
                }
                y = tmp;
            }

            if (z < 0)
            {
                int tmp = 0;
                int dir = 1;

                while (z < 0)
                {
                    tmp += dir;
                    if (tmp == height - 1 || tmp == 0)
                        dir *= -1;
                    z++;
                }
                z = tmp;
            }

            return data[getPos(x, y, z)];
        }

        public void set(float value, int x, int y, int z)
        {
            data[getPos(x, y, z)] = value;
        }

        public FloatArray2D getXPlane(int x)
        {
            FloatArray2D plane = new FloatArray2D(height, depth);

            for (int y = 0; y < height; y++)
                for (int z = 0; z < depth; z++)
                    plane.set(this.get(x, y, z), y, z);

            return plane;
        }

        public float[][] getXPlane_float(int x)
        {
            float[][] plane = new float[height][depth];

            for (int y = 0; y < height; y++)
                for (int z = 0; z < depth; z++)
                    plane[y][z] = this.get(x, y, z);

            return plane;
        }

        public FloatArray2D getYPlane(int y)
        {
            FloatArray2D plane = new FloatArray2D(width, depth);

            for (int x = 0; x < width; x++)
                for (int z = 0; z < depth; z++)
                    plane.set(this.get(x, y, z), x, z);

            return plane;
        }

        public float[][] getYPlane_float(int y)
        {
            float[][] plane = new float[width][depth];

            for (int x = 0; x < width; x++)
                for (int z = 0; z < depth; z++)
                    plane[x][z] = this.get(x, y, z);

            return plane;
        }

        public FloatArray2D getZPlane(int z)
        {
            FloatArray2D plane = new FloatArray2D(width, height);

            for (int x = 0; x < width; x++)
                for (int y = 0; y < height; y++)
                    plane.set(this.get(x, y, z), x, y);

            return plane;
        }

        public float[][] getZPlane_float(int z)
        {
            float[][] plane = new float[width][height];

            for (int x = 0; x < width; x++)
                for (int y = 0; y < height; y++)
                    plane[x][y] = this.get(x, y, z);

            return plane;
        }

        public void setXPlane(FloatArray2D plane, int x)
        {
            for (int y = 0; y < height; y++)
                for (int z = 0; z < depth; z++)
                    this.set(plane.get(y, z), x, y, z);
        }

        public void setXPlane(float[][] plane, int x)
        {
            for (int y = 0; y < height; y++)
                for (int z = 0; z < depth; z++)
                    this.set(plane[y][z], x, y, z);
        }

        public void setYPlane(FloatArray2D plane, int y)
        {
            for (int x = 0; x < width; x++)
                for (int z = 0; z < depth; z++)
                    this.set(plane.get(x, z), x, y, z);
        }

        public void setYPlane(float[][] plane, int y)
        {
            for (int x = 0; x < width; x++)
                for (int z = 0; z < depth; z++)
                    this.set(plane[x][z], x, y, z);
        }

        public void setZPlane(FloatArray2D plane, int z)
        {
            for (int x = 0; x < width; x++)
                for (int y = 0; y < height; y++)
                    this.set(plane.get(x, y), x, y, z);
        }

        public void setZPlane(float[][] plane, int z)
        {
            for (int x = 0; x < width; x++)
                for (int y = 0; y < height; y++)
                    this.set(plane[x][y], x, y, z);
        }

    }
}



