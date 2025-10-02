import { Text } from "@/components/ui/text";

export const Home = () => {
  return (
    <div className="w-11/12 max-w-5xl flex flex-col items-stretch gap-2">
      <Text size="t1" className="mt-48">
        AI for Brain Tumor Diagnosis
      </Text>
      <Text size="h2">CS4641 Fall 2025 Team 24</Text>
      <Text size="h3">
        Ritika Calyanakoti, Aaron Fan, Hannah Hsiao, Aaron Hung, Sagarika
        Srinivasan
      </Text>
      <Text size="h2" id="overview" className="mt-8">
        Overview
      </Text>
      <div className="grid grid-cols-1 md:grid-cols-[auto_400px] gap-12">
        <Text>
          Brain tumors are among the most difficult diseases to diagnose
          accurately, as MRI scans require careful interpretation and subtle
          differences can be overlooked. Our project explores the use of
          computer vision and machine learning to improve the accuracy and speed
          of brain tumor detection. Using a Kaggle dataset of brain MRI scans,
          we are building a multi-class classification model capable of
          identifying gliomas, meningiomas, pituitary tumors, or healthy cases.
          By testing models such as VGGNet, ResNet, and Vision Transformers, and
          applying preprocessing techniques like normalization and augmentation,
          our goal is to reduce diagnostic errors, support overworked
          clinicians, and ultimately improve patient outcomes through more
          consistent and reliable tumor classification.
        </Text>
        <img
          src="/pituitary1002.jpg"
          alt="MRI scan of human head with pituitary tumor"
          className="w-full h-full rounded-lg object-contain shadow-sm"
        />
      </div>
    </div>
  );
};
