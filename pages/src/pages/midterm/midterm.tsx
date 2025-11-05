import { Divider } from "@/components/ui/divider";
import { Image } from "@/components/ui/image";
import { Text } from "@/components/ui/text";
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table";
import pituitaryImage from "@/assets/pituitary1002.jpg";
import raw_data from "@/assets/raw_data_summary.png";
import hist_cluster_4 from "@/assets/hist_cluster_4.png";
import hist_example_4 from "@/assets/hist_example_4.png";
import rings_cluster_4 from "@/assets/rings_cluster_4.png";
import rings_example_4 from "@/assets/rings_example_4.png";
import vgg16_eval from "@/assets/vgg16_eval.png";
import vgg16_loss from "@/assets/vgg16_loss.png";
import vgg16_pretrained_eval from "@/assets/vgg16_pretrained_eval.png";
import vgg16_pretrained_loss from "@/assets/vgg16_pretrained_loss.png";
import gantt from "@/assets/GanttChart.xlsx.png";

export const Midterm = () => {
  const references = [
    `A. Iannessi, H. Beaumont, C. Aguillera, F. Nicol, and A.-S.
          Bertrand, “The ins and outs of errors in oncology imaging: the DAC
          framework for radiologists,” Frontiers in Oncology, vol. 14, Oct.
          2024, doi: https://doi.org/10.3389/fonc.2024.1402838.`,
    `D. N. Louis et al., “The 2021 WHO Classification of Tumors of the
          Central Nervous System: a summary,” Neuro-Oncology, vol. 23, no. 8,
          Jun. 2021, doi: https://doi.org/10.1093/neuonc/noab106.`,
    `M. Martucci et al., “Magnetic Resonance Imaging of Primary Adult
          Brain Tumors: State of the Art and Future Perspectives,”
          Biomedicines, vol. 11, no. 2, p. 364, Feb. 2023, doi:
          https://doi.org/10.3390/biomedicines11020364.`,
    `N. Melarkode, K. Srinivasan, S. M. Qaisar, and P. Plawiak,
          “AI-Powered Diagnosis of Skin Cancer: A Contemporary Review, Open
          Challenges and Future Research Directions,” Cancers, vol. 15, no. 4,
          p. 1183, Jan. 2023, doi: https://doi.org/10.3390/cancers15041183.`,
    `S. R. Chandana, S. Movva, M. Arora, and T. Singh, “Primary Brain
          Tumors in Adults,” American Family Physician, vol. 77, no. 10, pp.
          1423–1430, May 2008, Available:
          https://www.aafp.org/pubs/afp/issues/2008/0515/p1423.html`,
    "F. Mölder et al., “Sustainable Data Analysis with snakemake,” F1000Research, vol. 10, p. 33, Jan. 2021. doi:10.12688/f1000research.29032.1",
    "F. Gaillard, Y. Baba, and D. Bell, “MRI sequences (overview),” Radiopaedia.org, Jun. 2015. doi:10.53347/rid-37346 ",
  ];

  return (
    <div className="flex flex-col w-11/12 max-w-5xl items-stretch">
      <Text size="t2" className="mb-8">
        Midterm
      </Text>
      <Text size="h2" id="introduction">
        Introduction
      </Text>
      <Divider />
      <div className="clear-both">
        <Image
          src={pituitaryImage}
          alt="MRI scan of human head with pituitary tumor"
          className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
        />
        <Text>
          Brain tumors are some of the most aggressive diseases, requiring early
          and accurate diagnosis for effective treatment. MRI scans are the
          standard test for brain tumor diagnosis, simultaneously distinguishing
          between various types of cancers and limiting patient exposure to
          radiation [3]. MRIs also measure treatment response and are useful for
          oncologists to decide when to stop treatments or if other methods
          should be incorporated [3]. That said, MRI evaluations by clinicians
          are not fully accurate. Like all decision making processes, scan
          reading is subjective and prone to error [1]. Such situations
          necessitate healthcare professionals to seek forms of diagnosis with
          higher accuracy. The complexity of brain tumors further complicates
          diagnosis [2]. The tumor imaging space is extremely complicated with
          clinicians relying on MRIs to differentiate among numerous tumor
          classifications [5]. These cancers spread to different parts of the
          brain affecting the cerebrospinal fluid, soft tissue, and early and
          mid brain structures [5]. Dependence on MRI scan interpretation
          introduces potential error and creates a need for computational
          methods that can consistently and accurately diagnose scans. Recently,
          artificial intelligence and machine learning approaches have emerged
          as promising tools in medical research and diagnoses, becoming
          increasingly common with detecting skin cancer, facilitating drug
          development, and interpreting sensor data [4]. This project provides
          an opportunity to enhance diagnostic accuracy and overall clinical
          decision-making.
        </Text>
      </div>
      <Text size="h2" id="problem">
        Problem and Motivation
      </Text>
      <Divider />
      <Text>
        We are attempting to solve medical imaging diagnosis issues within the
        field of healthcare. Currently, doctors must manually scan images of
        their patients which is both a time consuming task and prone to human
        error, potentially leading to errors in life altering medical decisions.
        Early and accurate brain tumor detection can significantly alter
        treatment options and improve patients’ prognoses. To offload some of
        the burden placed on healthcare professionals and increase the accuracy
        of diagnoses, we plan to build a machine learning model that can perform
        multi-class classification of brain MRIs, will be able to analyze a
        large dataset of images quickly and generate consistent, accurate, and
        easy-to-miss predictions for glioma, meningioma, pituitary, or healthy
        cases.
      </Text>
      <Text size="h2" id="methods">
        Methods
      </Text>
      <Divider />
      <Text size="h3" id="dataset">
        Dataset
      </Text>
      <Text>
        For our dataset, we are using the{" "}
        <a
          href="https://www.kaggle.com/datasets/ishans24/brain-tumor-dataset?select=glioma"
          className="underline hover:opacity-50 transition-all"
          target="_blank"
        >
          Brain Tumor Dataset
        </a>{" "}
        from Kaggle.
      </Text>
      <Text size="h3" id="preprocessing">
        Preprocessing
      </Text>
      <ul className="list-disc list-outside ml-8 space-y-2 pb-2">
        <li>
          <Text>
            <b>Grayscale</b>: The dataset we are using is assembled from several
            different datasets, each of which has different coloring. To make
            our dataset uniform across all samples, we can turn all images into
            grayscale to ensure the model is learning the patterns in the image
            and not the colors to classify the image.
          </Text>
        </li>
        <li>
          <Text>
            <b>Normalize</b>: Since different images can also have different
            pixel brightness ranges, we can also normalize all of our image
            pixel brightness to be between 0 and 1. This improves model
            stability by ensuring the model is resistant against variations in
            the data.
          </Text>
        </li>
        <li>
          <Text>
            <b>Image rescaling</b>: To make sure all of the images are the same
            dimensions, we scaled all images to follow the most common image
            size in our dataset, which is 512-by-512 pixels. The images were
            scaled so the smaller dimension is 512 px, and the other dimension
            is cropped. Some of the images were further scaled before training
            and running inference depending on the requirements of the model.
          </Text>
        </li>
        <li>
          <Text>
            <b>Data Augmentation</b>: We can rotate or flip our images,
            introduce noise, or blur/sharpen the image to reproduce common
            variations in image quality and allows the model to take those
            factors into account.
          </Text>
        </li>
      </ul>
      <Text>
        The data was also split into a training, cross-validation (CV), and
        testing set. The CV and testing sets were each 15% of the total data
        set. Each of 4 classes, after preprocessing, had a total of 3754 images,
        which means 563 images were used for CV and for testing. Each set is
        mutually exclusive, so none of the images overlapped.
      </Text>
      <Text size="h3" id="preprocessing">
        Machine Learning Algorithms
      </Text>
      <Text>
        We also have three potential machine learning models (one has been
        implemented so far) we can test to select the one that performs the
        best:
      </Text>
      <ul className="list-disc list-outside ml-8 space-y-2">
        <li>
          <Text>
            <b>VGGNet</b>: This is a simple CNN that only uses 3x3 convolutional
            filters but still produces high-quality results. This serves as our
            baseline model to compare the more complex models against to see how
            model complexity affects performance.
          </Text>
        </li>
        <li>
          <Text>
            <b>ResNet</b>: ResNet adds residual connections to the model,
            allowing the network to skip layers and alleviate the vanishing
            gradient problem. This enables the training of much deeper networks
            than traditional CNNs while maintaining high accuracy, and exploring
            whether deeper models can capture more complex features.
          </Text>
        </li>
        <li>
          <Text>
            <b>Swin Transformer</b>: This model leverages a hierarchical vision
            transformer architecture that divides an image into non-overlapping
            windows and computes self-attention within each window to capture
            contextual information and complex patterns.
          </Text>
        </li>
      </ul>
      <Text size="h2" id="results">
        Training Process
      </Text>
      <Divider />
      <Text>
        The training process for the models were all the same. During training,
        the data is shuffled and split into batches of 64 images to train the
        model on. The cross entropy loss function was used to evaluate the model
        and perform backpropagation, and the Adam optimizer with a learning rate
        of 0.0001 was used to train the models. After each epoch, the model was
        evaluated against the cross validation set as a predictor of the
        training accuracy. The models were trained until the loss didn’t
        decrease for at least 10 epochs, or up to 100 epochs. The CV dataset was
        used instead of the testing dataset to avoid confirmation bias.
        Afterwards, each model was evaluated against the testing dataset, and
        metrics for the confusion matrix, F-score, and AUC-ROC were calculated.
      </Text>
      <Text>
        The initial exploratory data analysis, data preprocessing, model
        training, and model evaluation were all performed on PACE ICE using a
        pipeline set up via Snakemake [6].
      </Text>
      <Text size="h2" id="results">
        Results and Discussion
      </Text>
      <Divider />
      <Text size="h3" id="metrics">
        Exploratory Data Analysis
      </Text>

      <Text>
        Before data preprocessing and data analysis, the data analysis was
        performed to get an idea of the contents of the data and any potential
        failure points for the models. Figure 1 is the results of this initial
        data analysis.
      </Text>

      <Image
        src={raw_data}
        alt="Figure 1: Summary of image dimensions and image contrasts for each class across the dataset. Top: 8 most common image dimensions by image class. Bottom: violin plot of image contrasts by image class."
        className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
      />

      <Text>
        From our data analysis, 2 potential issues were identified. First, the
        dimensions for all images are relatively similar, except for the
        no_tumor class, which had a wide variety of image dimensions and none
        that were the same as the rest of the classes, and that the number of
        images in each class are far from evenly distributed. This discrepancy
        meant data augmentation was needed to prevent one specific class from
        being more common than the rest, resulting in poor model performance
        resulting from optimization of just the most common class. Second, the
        image contrasts were all generally pretty similar, except again for the
        no tumor class. The no tumor class had a significantly higher contrast
        than the rest of the classes, despite having no tumor. This is likely
        due to the issue of healthy samples being rare in medical data: if
        patients are healthy, why bother collecting their data? This could pose
        a potential issue with our model if it ends up learning that higher
        contrast indicates no tumor, which is not a result we want to see.
        Another interesting point from the second figure is the bimodal
        distribution of the contrasts of the other 3 classes. Upon further
        inspection, this seemed to be due to the presence of different MRI
        sequences in the data. Depending on the sequence that radio frequencies
        and pulses are sent by the MRI machine, different types of tissue within
        the brain will be highlighted in different colors [7]. This means that
        different images within the same class can look significantly different,
        which can potentially lead to poor results. This is likely what caused
        the bimodal distribution between the image contrasts of each class.
        Although the images cannot be converted from one MRI sequence to
        another, one potential solution was to perform simple clustering to
        identify the sequences for each image, and pass that information to the
        model to give it context.
      </Text>

      <Text size="h3" id="metrics">
        MRI Sequence Clustering
      </Text>
      <Text>
        To perform clustering for the MRI sequences, we performed K-means
        clustering with two different types of engineered features. First, since
        the contrast was what was causing the split distribution, we thought we
        might be able to cluster the pixels based on their pixel intensity
        histograms. The second method was using the knowledge that since the
        brain is roughly circular and symmetrical, we can use average pixel
        brightnesses at different radii around the image to cluster the images.
        We performed clustering on both methods from 1 to 18 clusters,
        calculated the distortion for each clustering, and used the elbow method
        to identify a suitable number of clusters. The results can be seen in
        figure 2 and 3.
      </Text>

      <Image
        src={hist_cluster_4}
        alt="Figure 2: Distortion vs. number of clusters for brightness histogram clustering"
        className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
      />
      <Image
        src={rings_cluster_4}
        alt="Figure 3: Distortion vs. number of clusters for clustering based on average radial pixel brightness"
        className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
      />

      <Text>
        For both methods, a good number of clusters based on the elbow method is
        4 clusters. Although it might seem from figure 3 that clustering based
        on average radial pixel brightness results in lower distortion, the
        difference in scale is largely due to the brightness histograms having a
        higher number of dimensions than the average radial pixel brightness. To
        be able to precisely evaluate whether either clustering method produces
        useful information, an example for each cluster from each class was also
        plotted in figures 4 and 5.
      </Text>

      <Image
        src={hist_example_4}
        alt="Figure 4: Example images for clustering based on pixel brightness histograms"
        className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
      />
      <Image
        src={rings_example_4}
        alt="Figure 5: Example images for clustering based on average radial pixel brightness"
        className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
      />

      <Text>
        From figures 4 and 5, it is clear that neither method was able to truly
        differentiate between different MRI sequences. From these images, it
        seems to be because the images may be dominated by a few MRI sequences,
        resulting in the clustering not being able to precisely capture the
        variations in pixel brightness for each sequence. Because of this, we
        decided it was not worth pursuing this avenue for feature engineering
        and just to use the data that is already present.
      </Text>

      <Image
        src={vgg16_loss}
        alt=""
        className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
      />

      <Text>
        The above graph is classic training behavior. The loss starts at a high
        value, around 31, and then drops rapidly as the network learns.
        represents the Training loss vs. Epoch for the VGGNet which is our
        baseline model. As we can see the number of epochs does not drastically
        improve after 5 epochs and keeps oscillating within the 4-7 loss range.
        After epoch 15 we can say the model converged and that additional
        training does not significantly improve performance.
      </Text>

      <Image
        src={vgg16_pretrained_loss}
        alt=""
        className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
      />

      <Text>
        The above graph has a slightly unusual curve. The model starts at a
        slightly lower loss value, 29, which is expected given that the
        pre-trained data is less random. That said, there is a spike around
        epoch 12 and then a sharp decline around 14. It could be an example of
        gradient exploding which is where the weight updates become so large
        that the network is destabilized, causing the loss to drastically
        increase.
      </Text>

      <Image
        src={vgg16_eval}
        alt=""
        className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
      />

      <Text>
        The Confusion Matrix allows us to evaluate how well the model worked.
        The values in the boxes tell us how many models were matched to a
        specific label. In this case, we are looking at the vgg16 model. By
        looking at the dark blue diagonals, we see that the majority of the
        images were mapped to the correct label. The F-score is 0.96130,
        highlighting that there is both high precision and high recall across
        specific thresholds. This means that our model has few false positives
        and few false negatives. Additionally, the AUC being 0.99316 highlights
        that our model can properly distinguish between positive and negative
        classes between all possible thresholds, meaning that our model has been
        able to properly learn how to distinguish between the images.
      </Text>

      <Image
        src={vgg16_pretrained_eval}
        alt=""
        className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
      />

      <Text>
        Here, we are running the vgg16_pretrained model. Again, the diagonal of
        the confusion matrix holds a majority of the values, highlighting that
        once again, the majority of the images were mapped to the correct label.
        The F-score is 0.95979, highlighting that there is both high precision
        and high recall across specific thresholds. This means that our model
        has few false positives and few false negatives. Additionally, the AUC
        being 0.99623 highlights that our model can properly distinguish between
        positive and negative classes between all possible thresholds, meaning
        that our model has been able to properly learn how to distinguish
        between the images.
      </Text>

      <Text size="h3" id="next-steps">
        Next Steps
      </Text>
      <Text>
        Our next steps are to complete training and evaluation of the remaining
        2 models, and to compare the results of each model. Since our first
        naive model already seems to perform quite well, these remaining, more
        complex models are also expected to perform at a similar level, so we
        might also attempt using even simpler models to figure out what the
        simplest model is that can still effectively differentiate between the
        different types of brain tumor.
      </Text>

      <Text size="h2" id="responsibilities">
        Responsibilities
      </Text>
      <Divider />
      <Text size="h3" id="gantt-chart">
        Gantt Chart
      </Text>
      <img
        src={gantt}
        alt="Gantt chart of planned tasks duration and deadlines"
        className="m-4"
      />
      <Text size="h3" id="contributions">
        Contribution Table
      </Text>
      <Table className="mt-4">
        <TableHeader>
          <TableRow>
            <TableHead>Name</TableHead>
            <TableHead>Proposal Contributions</TableHead>
          </TableRow>
        </TableHeader>
        <TableBody>
          <TableRow>
            <TableCell>Aaron Hung</TableCell>
            <TableCell>Report, Model Training</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Aaron Fan</TableCell>
            <TableCell>Report, Data Sourcing + Cleaning</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Ritika Calyanakoti</TableCell>
            <TableCell>Report, Data Sourcing + Cleaning</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Hannah Hsiao</TableCell>
            <TableCell>Report, Model Training</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Sagarika Srinivasan</TableCell>
            <TableCell>Report, Model Coding</TableCell>
          </TableRow>
        </TableBody>
      </Table>
      <Text size="h2" id="references">
        References
      </Text>
      <Divider />
      <ol>
        {references.map((ref, i) => (
          <li className="flex flex-row gap-1">
            <Text>[{i + 1}]</Text>
            <Text>{ref}</Text>
          </li>
        ))}
      </ol>
    </div>
  );
};
