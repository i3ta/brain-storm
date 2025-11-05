import { Divider } from "@/components/ui/divider";
import { Image } from "@/components/ui/image";
import { Text } from "@/components/ui/text";
import Latex from "react-latex-next";
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table";
import pituitaryImage from "@/assets/pituitary1002.jpg";
import gantt from "@/assets/GanttChart.xlsx.png";

export const Proposal = () => {
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
  ];

  return (
    <div className="flex flex-col w-11/12 max-w-5xl items-stretch">
      <Text size="t2" className="mb-8">
        Proposal
      </Text>
      <Text size="h2" id="introduction">
        Introduction
      </Text>
      <Divider />
      <div className="clear-both">
        <Image
          float
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
            <b>Data Augmentation</b>: We can rotate or flip our images,
            introduce noise, or blur/sharpen the image to reproduce common
            variations in image quality and allows the model to take those
            factors into account.
          </Text>
        </li>
      </ul>
      <Text size="h3" id="preprocessing">
        Machine Learning Algorithms
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
        Results
      </Text>
      <Divider />
      <Text size="h3" id="metrics">
        Quantitative Metrics
      </Text>
      <ul className="list-disc list-outside ml-8 space-y-2 pb-2">
        <li className="clear-both">
          <Image
            float
            src="https://media.geeksforgeeks.org/wp-content/uploads/20240708132251/confusion-Matrix.PNG"
            alt="Image of confusion matrix of 3 classes. Source: GeeksForGeeks"
            className="max-w-sm ml-4 mb-4"
          />
          <Text>
            <b>Confusion Matrix</b>: This is a table of true positive, true
            negative, false positive, and false negative rates for each class.
            In a medical context where we might prefer false positives over
            false negatives, a confusion matrix helps us evaluate the accuracy
            of our model.
          </Text>
        </li>
        <li>
          <div className="flex flex-col items-stretch gap-2">
            <Text>
              <b>F1 Score</b>: The F1 score is a metric balanced between
              precision and recall into one value as follows:
            </Text>
            <div className="flex justify-center">
              <Latex>
                {`$F_1 = \\frac{2 \\times \\text{Precision} \\times \\text{Recall}}{\\text{Precision} + \\text{Recall}}$`}
              </Latex>
            </div>
            <Text>
              This allows us to compare the performance of each model directly
              without worrying about imbalances between false positives and
              false negatives, making it especially useful with datasets where
              classifications may be underrepresented.
            </Text>
          </div>
        </li>
        <li>
          <Text>
            <b>AUC-ROC Score</b>: This score, well-established in medical
            imaging, is calculated by integrating the area under the receiver
            operating characteristic curve, which is a plot of the true positive
            rate against the true negative rate. This allows us to evaluate
            different hyperparameters to find the best for tumor identification.
          </Text>
        </li>
      </ul>
      <Text size="h3" id="goals">
        Project Goals
      </Text>
      <Text>
        Ultimately, our goal is to increase accuracy when diagnosing patients’
        brain scans for tumors. We want to streamline medical processes so that
        patient care and treatment is the top priority, rather than manually
        examining images and consulting other professionals for second opinions.
      </Text>
      <Text size="h3" id="expectations">
        Expected Results
      </Text>
      <Text>
        Our model is expected to classify an image of a brain scan and classify
        it as glioma, meningioma, pituitary, or healthy.
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
            <TableCell>Methods, Slides</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Aaron Fan</TableCell>
            <TableCell>Slides, Results</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Ritika Calyanakoti</TableCell>
            <TableCell>Problem Definition, Slides</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Hannah Hsiao</TableCell>
            <TableCell>Data Overview, Gantt Chart</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Sagarika Srinivasan</TableCell>
            <TableCell>Introduction and Background</TableCell>
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
