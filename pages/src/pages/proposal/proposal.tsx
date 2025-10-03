import { Divider } from "@/components/ui/divider";
import { Image } from "@/components/ui/image";
import { Text } from "@/components/ui/text";
import Latex from "react-latex-next";
import pituitaryImage from "@/assets/pituitary1002.jpg";

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
          src={pituitaryImage}
          alt="MRI scan of human head with pituitary tumor"
          className="h-full max-h-96 rounded-lg object-contain shadow-sm m-4"
        />
        <Text>
          Brain tumors are some of the most aggressive diseases, requiring early
          and accurate diagnosis for effective treatment. MRI scans are the
          standard test for brain tumor diagnosis simultaneously distinguishing
          between various types of cancers and limiting patient exposure to
          radiation [3]. MRIs also measure treatment response and are useful for
          oncologists to decide when to stop treatments, if other methods should
          be incorporated, etc. [3]. That being said, MRI evaluations by
          clinicians are not fully accurate. Oncologists frequently make errors
          in interpreting cancer imaging such as misreading scans or having
          their own cognitive biases [1]. Such situations necessitate healthcare
          professionals to seek forms of diagnosis with higher accuracy. The
          complexity of brain tumors further complicates diagnosis. The tumor
          imaging space is extremely complicated with clinicians relying on MRs
          to differentiate among numerous tumor classifications. Examples
          include Gliomas, Meningiomas, Lymphomas, Pineal tumors, etc [5]. These
          cancers spread to different parts of the brain with some affecting the
          cerebrospinal fluid, others the soft tissue, and even more affecting
          the early and mid brain structures [5]. With so many tumor types and
          affected regions, it is important to have a reference or baseline
          tumor aggressiveness/intensity scale. There are existing solutions for
          this. For example, the World Health Organization's tumor grading
          system is the widely accepted standard which ranks the central nervous
          system tumors from Grade I to Grade IV (least to most malignant) [2].
          While the WHO grading system is a standardized framework, its
          application doesn’t negate reliance on clinicians’s interpretation of
          MRI scans. This dependence introduces potential error and creates a
          need for computational methods that can improve consistency and
          accuracy in diagnosis. As of the past decades, artificial intelligence
          and machine learning approaches have emerged as promising tools in
          medical imaging. Machine Learning applications in healthcare are now
          increasingly common with detecting rates of skin cancer, drug
          development, human-monitoring sensors [4]. This project provides an
          exciting opportunity to enhance diagnostic accuracy and overall
          clinical decision-making.
        </Text>
      </div>
      <Text size="h2" id="problem">
        Problem and Motivation
      </Text>
      <Divider />
      <Text>
        We are attempting to solve medical imaging diagnosis issues within the
        field of healthcare. Currently, doctors must manually scan images of
        their patients which is both a time consuming task, and one prone to
        human error. Diagnoses might vary from doctor to doctor since
        interpretations of brain scans vary. This can lead to potential errors
        in life altering medical decisions, as many cases might be overlooked,
        such early-stage tumors, or falsely diagnosed. Early and accurate brain
        tumor detection can significantly alter treatment options and improve
        patients’ prognoses, creating a critical need for a more accurate
        methodology. In today’s healthcare industry, there is also a shortage of
        doctors, leading hospitals to be overworked and understaffed. To offload
        some of the burden placed on healthcare professionals and increase the
        accuracy of diagnoses, we plan to build a machine learning model. This
        model, used for multi-class classification of brain MRIs, will be able
        to analyze a large dataset of images quickly and generate consistent,
        accurate, and easy-to-miss predictions for glioma, meningioma,
        pituitary, or healthy cases.
      </Text>
      <Text size="h2" id="methods">
        Methods
      </Text>
      <Divider />
      <Text>
        To find a model that performs best for this task, we plan on employing 3
        preprocessing methods and evaluating 3 separate machine learning models.
      </Text>
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
      <Text>To preprocess our data, we have three potential approaches:</Text>
      <ul className="list-disc list-outside ml-8 space-y-2 pb-2">
        <li>
          <Text>
            <b>Grayscale</b>: The dataset we are using is assembled from several
            different datasets, each of which has slightly different coloring.
            To make our dataset uniform across all samples, we can turn all
            images into grayscale to ensure the model is learning the patterns
            in the image, not just the colors of the image to determine whether
            there is a tumor or not.
          </Text>
        </li>
        <li>
          <Text>
            <b>Normalize</b>: Since different images can also have different
            pixel brightness ranges, especially after making all images
            grayscale, we can also normalize all of our image pixel brightness
            to be between 0 and 1. This also helps improve model stability by
            ensuring the model is less affected by variations in the scale of
            the data.
          </Text>
        </li>
        <li>
          <Text>
            <b>Data Augmentation</b>: Although our dataset is already quite
            large, we can still perform data augmentation to help improve the
            model's robustness to different images. To do this, we can rotate or
            flip our images, introduce small amounts of noise, or blur/sharpen
            the image to reproduce common variations in image quality and allows
            the model to take those factors also into account.
          </Text>
        </li>
      </ul>
      <Text size="h3" id="preprocessing">
        Machine Learning Algorithms
      </Text>
      <Text>
        We also have three potential machine learning models we can test to
        select the one that performs the best:
      </Text>
      <ul className="list-disc list-outside ml-8 space-y-2">
        <li>
          <Text>
            <b>VGGNet</b>: This is a simple convolutional neural network (CNN)
            that only uses 3x3 convolutional filters but has been reported to
            still produce high-quality results. This serves as our baseline
            model to compare the other, more complex models against to see how
            much increasing the complexity of the model affects the performance,
            if at all.
          </Text>
        </li>
        <li>
          <Text>
            <b>ResNet</b>: ResNet introduces the concept of residual
            connections, which allow the network to skip layers and prevent the
            vanishing gradient problem. This architecture enables the training
            of much deeper networks compared to traditional CNNs while
            maintaining high accuracy. Using ResNet allows us to explore whether
            deeper models with residual learning can capture more complex
            features in the images compared to simpler networks like VGGNet.
          </Text>
        </li>
        <li>
          <Text>
            <b>Swin Transformer</b>: This model leverages a hierarchical vision
            transformer architecture that divides an image into non-overlapping
            windows and computes self-attention within each window. By doing so,
            it captures both local and global contextual information
            efficiently, often outperforming traditional CNNs on image
            recognition tasks while using fewer parameters. This makes it a
            strong candidate for capturing complex patterns in medical imaging.
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
      <Text>
        To evaluate our models, we plan on using the following metrics:
      </Text>
      <ul className="list-disc list-outside ml-8 space-y-2 pb-2">
        <li className="clear-both">
          <Image
            src="https://media.geeksforgeeks.org/wp-content/uploads/20240708132251/confusion-Matrix.PNG"
            alt="Image of confusion matrix of 3 classes. Source: GeeksForGeeks"
            className="max-w-sm ml-4 mb-4"
          />
          <Text>
            <b>Confusion Matrix</b>: This is a table of true positive, true
            negative, false positive, and false negative rates for each class.
            In a medical context where we might prefer false positives over
            false negatives, a confusion matrix would helps us evaluate the
            accuracy of our model in each case.
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
              This allows us to compare the performance of each model directly,
              and not have to worry about imbalances between false positives and
              false negatives, making it especially useful when dealing with
              datasets where certain classes (like tumor regions) may be
              underrepresented.
            </Text>
          </div>
        </li>
        <li>
          <Text>
            <b>AUC-ROC Score</b>: This is a score calculated by integrating the
            area under the receiver operating characteristic (ROC) curve, which
            is just a curve that plots the true positive rate against the true
            negative rate. This allows us to evaluate the performance of
            different variations of the same model with different
            hyperparameters together to find the model best suited for the
            specific situation. This metric is also useful since it has been
            used for medical imaging and has well-established thresholds that
            can be followed in determining the effectiveness of the model and
            whether it can be used for medical applications.
          </Text>
        </li>
      </ul>
      <Text size="h3" id="goals">
        Project Goals
      </Text>
      <Text>
        Our goal is to increase accuracy when diagnosing patients’ brain scans
        for tumors. This improves the ethical nature in medicine as patients
        will no longer be at risk of having a doctor misdiagnose them.
        Ultimately, we want to streamline medical processes so that patient care
        and treatment is the top priority, rather than manually examining images
        and consulting other professionals for second opinions.
      </Text>
      <Text size="h3" id="expectations">
        Expected Results
      </Text>
      <Text>
        Our model is expected to classify an image of a brain scan and classify
        it as glioma, meningioma, pituitary, or healthy.
      </Text>
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
