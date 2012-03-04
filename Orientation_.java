/* Orientation measurement plugin for ImageJ
 * 
 * This plugin uses statistics to calculate an orientation vector
 * for each label (segment) in the image. Stacks are considered 3D images.
 * 
 * Theory of operation:
 * Calculates the orientation of all labels (intensity values) in the image
 * by finding the eigenvector associated with the largest eigenvalue of
 * the covariance matrix of the x,y,z coordinates of each label.
 * Jama (Java matrix package, public domain) is used for extracting the 
 * eigenvalues and eigenvectors from the covariance matrix.
 * The calculation of the covariance matrix uses the fact that 
 * S_ij = E{Xi-E[Xi][Xj-E[Xj]} = E[Xi*Xj] -E[Xi]*E[Xj]
 * A map is used for storing data about each label in order to 
 * support sparse labels.
 * The image is read only once. The complexity of the algorithm 
 * is O(N), where N is the number of pixels in the image.
 * 
 * Author: Per Christian Henden, perchrh [at] pvv.org
 * License: public domain. Free for all use.
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageWindow;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.text.DecimalFormat;
import java.util.HashMap;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class Orientation_ implements PlugInFilter {

	private class CovarianceData {
		private long N = 0;

		private long sum_x = 0;

		private long sum_y = 0;

		private long sum_z = 0;

		private long sum_xx = 0;

		private long sum_zz = 0;

		private long sum_yy = 0;

		private long sum_xy = 0;

		private long sum_xz = 0;

		private long sum_yz = 0;

		private Matrix covariance = null;

		public void update(int x, int y, int z) {
			N++;
			sum_x += x;
			sum_y += y;
			sum_z += z;

			sum_xx += x * x;
			sum_yy += y * y;
			sum_xy += x * y;

			sum_zz += z * z;
			sum_yz += y * z;
			sum_xz += x * z;
		}

		public Matrix getMatrix(int dimension) {

			double x_bar = sum_x / N; // E[X]
			double y_bar = sum_y / N; // E[Y]
			double z_bar = sum_z / N; // E[Z]

			double cov_xx = sum_xx / N - x_bar * x_bar;
			double cov_yy = sum_yy / N - y_bar * y_bar;
			double cov_zz = sum_zz / N - z_bar * z_bar;

			double cov_xy = sum_xy / N - y_bar * x_bar;
			double cov_xz = sum_xz / N - x_bar * z_bar;
			double cov_zy = sum_yz / N - z_bar * y_bar;

			if (dimension == 2) {
				double[] tmp = { cov_xx, cov_xy, cov_xy, cov_yy };
				covariance = new Matrix(tmp, dimension);
			} else {
				double[] tmp = { cov_xx, cov_xy, cov_xz, cov_xy, cov_yy,
						cov_zy, cov_xz, cov_zy, cov_zz };
				covariance = new Matrix(tmp, dimension);
			}
			return covariance;
		}

		public double[] getCentroid(){
			double centerX = (double)(sum_x / N );
			double centerY = (double)(sum_y / N );

			return new double[] { centerX, centerY };
		}

	}

	private ImagePlus imp;
	private boolean doPlot;
	private boolean noGo;


	public void run(final ImageProcessor ip) {

		if (noGo){
			return;
		}

		final ImageStack stack = imp.getStack();
		final int depth = imp.getStackSize();
		final int height = imp.getHeight();
		final int width = imp.getWidth();
		ImageProcessor sip = null;
		HashMap<Integer, CovarianceData> labelData = new HashMap<Integer, CovarianceData>();

		int dimension = 2;
		if (depth > 1) {
			dimension = 3;
		}

		for (int z = 1; z <= depth; z++) {
			IJ.showProgress(z, depth);
			sip = stack.getProcessor(z);
			for (int x = 0; x < width; x++) {
				for (int y = 0; y < height; y++) {
					int label = sip.get(x, y);

					if (!labelData.containsKey(label)) {
						labelData.put(label, new CovarianceData());
					}

					labelData.get(label).update(x, y, z);
				}
			}
		}

		IJ.write("Results of orientation analysis of image " + imp.getTitle() + ":");
		for (int label : labelData.keySet()) {
			Matrix covariance = labelData.get(label).getMatrix(dimension);
			double[] orientation = getOrientation(covariance);
			IJ.write("Class " + label + " has orientation : " 
					+ printMatrix(orientation));
		}

		if (doPlot){
			IJ.write("Plotting orientation vectors to new image...");
			drawVectors(labelData);
		}
	}

	private double[] getOrientation(Matrix covariance) {
		int dimension = covariance.getColumnDimension();
		// Get the largest eigenvalue and its associated eigenvector
		EigenvalueDecomposition ed = covariance.eig();
		double[] orientation = new double[dimension];
		for (int col = 0; col < dimension; col++) {
			// The largest eigenvalue's eigenvector is found at the final row
			orientation[col] = ed.getV().get(col, dimension - 1);
		}
		// Normalize the vector
		return normalizeVector(orientation);
	}



	private double[] normalizeVector(double[] input) {
		double sum = 0;
		final int length = input.length;
		for (int i = 0; i < length; i++) {
			sum += input[i] * input[i];
		}
		sum = Math.sqrt(sum);
		double[] output = new double[length];
		for (int i = 0; i < length; i++) {
			output[i] = input[i] / sum;
		}
		return output;
	}

	private String printMatrix(double[] input) {
		double[][] wrap = new double[1][];
		wrap[0] = input;
		return printMatrix(wrap);
	}

	private String printMatrix(double[][] input) {
		final DecimalFormat df5 = new DecimalFormat("#0.00000");
		StringBuffer out = new StringBuffer();
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input[i].length; j++) {
				out.append(df5.format(input[i][j]));
				if (j != input[i].length - 1) {
					out.append(", ");
				}
			}
			out.append("\n");
		}
		return new String(out);
	}


	/* ImageJ specific functions */

	private void drawVectors(HashMap<Integer, CovarianceData> labelData){

		final int centroidWidth = 6;
		final int strokeWidth = 1;
		final int vectorLength = 75;

		try{
			ImagePlus thePlot = new ImagePlus();
			thePlot.setProcessor("Orientation vector display", imp.getProcessor().convertToRGB());
			Image theImage = thePlot.getImage();

			Graphics2D g2d = (Graphics2D)theImage.getGraphics();
			g2d.setPaint(Color.magenta);
			g2d.setStroke(new BasicStroke(strokeWidth));
			for (int label : labelData.keySet()) {
				double[] centroid = labelData.get(label).getCentroid();
				Matrix covariance = labelData.get(label).getMatrix(2);
				double[] orientation = getOrientation(covariance);
				Point2D startPoint = new Point2D.Double(centroid[0], centroid[1]);
				Point2D endPoint = new Point2D.Double(startPoint.getX() + vectorLength * orientation[0], startPoint.getY() + vectorLength * orientation[1]);

				Ellipse2D theCentroid = new Ellipse2D.Double(centroid[0]-centroidWidth/2, centroid[1]-centroidWidth/2, centroidWidth, centroidWidth);
				g2d.fill(theCentroid);
				g2d.draw(theCentroid);

				Line2D theVector = new Line2D.Double(startPoint, endPoint);
				g2d.fill(theVector);
				g2d.draw(theVector);

			}

			ImageWindow f = new ImageWindow(thePlot);
			f.setEnabled(true);
		} catch (Exception e){

			IJ.write("There was a problem visualizing the vectors. The problem was:");
			IJ.write(e.toString());
		}

		return;
	}


	public int setup(String arg, final ImagePlus imp) {
		this.imp = imp;
		noGo = false;
		if (arg.equals("about")) {
			showAbout();
			return DONE;
		}

		if (imp.getStackSize() == 1){
			doPlot = true; // Plot results in additional image when 2D
		}

		if (!noGo) {
			IJ.showStatus("Measuring orientation..");
			IJ.write("Measuring orientation..");
		}

		return DOES_16 + DOES_8G + NO_CHANGES;
	}


	private void showAbout() {
		IJ.showMessage("Calculates the orientation of all labels (intensity values) in the image\n"
				+ "by finding the eigenvector associated with the largest eigenvalue of\n"
				+ "the covariance matrix of the x,y,z coordinates of each label.");
	}
}
