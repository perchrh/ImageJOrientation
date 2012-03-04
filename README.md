ImageJ plugins::Orientation measurement
====================

This plugin uses statistics to calculate an orientation vector
for each label (segment) in the image. 
This vector describes the direction (orientation) of the segment.
Stacks are considered 3D images. 
A plot of the orientation vectors is shown.

[Download ready to run plugin](http://www.pvv.org/~perchrh/imagej/Orientation-1.0.1_.jar)

Example output
---------------------


Results of orientation analysis of image fibers149-01.jpg:

....

Class 44 has orientation : -0.10987, 0.99395

Class 45 has orientation : -0.37344, 0.92765

Class 46 has orientation : -0.14266, 0.98977

Class 47 has orientation : 0.95784, -0.28731

....


How to install
---------------------

Download the plugin and copy it to your ImageJ plugin directory.
Java 1.5 or later is required to use the plugin.

How to compile it yourself
---------------------

- Using Java 1.5 or later, with Jama on the classpath, compile \*.java
- Build a jar file Using Java 1.5 or later including \*.class and plugin.properties on the root directory and a Jama/ directory containing the Jama .class files


Theory of operation
---------------------

Calculates the orientation of all labels (intensity values) in the image
by finding the eigenvector associated with the largest eigenvalue of
the covariance matrix of the x,y,z coordinates of each label.
For the mathematical details, see [this paper](http://www.pvv.org/~perchrh/papers/mastersthesisHendenBacheWiig.pdf), pages 31-32 (23-24).

[Jama](http://math.nist.gov/javanumerics/jama/), the Java matrix package, public domain, is used for extracting the eigenvalues and eigenvectors from the covariance matrix.  
The vectors are normalized before printed. 

The coordinate system of the image is used. That means (x=0,y=0) is at the topmost left corner of the screen. X-coordinates increase towards the right and y-coordinates increase downwards. Z-coordinates increase into the screen.

The image is read only once. The complexity of the algorithm is O(N), where N is the number of pixels in the image.

License: public domain. Free for all use.
