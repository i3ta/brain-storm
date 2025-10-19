## Exporatory Data Analysis Notes

Started out plotting the image dimensions and the contrast of the images for all
image classes:

![Raw Data Summary](https://github.gatech.edu/chung62/cs4641-team24/blob/a8c9ae21960c80335bdd2ac1aa86dc904f9a7f55/workflows/results/brain-tumor-dataset/summary/raw_data_summary.png)

From this, it was clear there were 2 things that needed to be addressed:

1. Image dimensions and image counts across classes
2. Image contrast

The second one is more concerning. The first one can be solved simply by adding
padding and using data augmentation. However, this figure isn't very descriptive
when it comes to the distribution of the image sizes, so a second plot was made
to highlight that.

![Raw image dimensions](https://github.gatech.edu/chung62/cs4641-team24/blob/a8c9ae21960c80335bdd2ac1aa86dc904f9a7f55/workflows/results/brain-tumor-dataset/summary/raw_image_dims.png)

From this it is clear that all of the classes have most images that are 512 x
512, except for the no tumor class that seems to have a variety of image
dimensions. That class also has fewer images in general.

From the figure, all 3 classes have similar root mean square contrasts, except
for the no tumor class. This is more concerning, because if the images are
different then the model may be learning the wrong information, leading to data
leakage. The model may be learning that higher contrast implies no tumor, which
is not what the model should be looking for. However, this may also be due to
the images having different MRI sequences (T1, T2, FLASK) that result in
different types of brain matter being highlighted different colors.

To deal with this, we tried using K-Means clustering to see if it is possible to
separate the different MRI sequences. To make this simpler, we calculated a
normalized histogram of all the pixel intensities, and performed K-Means
clustering for all values from 1 to 18 clusters (18 being chosen because we
thought with 3 main types of MRI sequences and 4 orientations (top, front, side,
back) we might end up with around 12 clusters, and 1.5 times that would allow us
to verify). We plotted the number of clusters against the distortion, and got
the following plot.

![Distortion chart with 4 clusters](https://github.gatech.edu/chung62/cs4641-team24/blob/a8c9ae21960c80335bdd2ac1aa86dc904f9a7f55/workflows/results/brain-tumor-dataset/hist/cluster_4.png)

Using the elbow method, we decided that 4 clusters is a good place to start.
From getting example images for each of the clusters, we got the following:

![Example images with 4 clusters](https://github.gatech.edu/chung62/cs4641-team24/blob/a8c9ae21960c80335bdd2ac1aa86dc904f9a7f55/workflows/results/brain-tumor-dataset/hist/example_4.png)

From this image, it is pretty clear that the clusters do not correspond to the
different types of MRI sequences, so we tested out more clusters. Below is the
example with 12 clusters:

![Example images with 12 clusters](https://github.gatech.edu/chung62/cs4641-team24/blob/a8c9ae21960c80335bdd2ac1aa86dc904f9a7f55/workflows/results/brain-tumor-dataset/hist/example_12.png)

This clustering still did not provide a good separation between MRI sequences.
