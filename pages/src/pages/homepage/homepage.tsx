import { Text } from "@/components/ui/text";
import pituitaryImage from "@/assets/pituitary1002.jpg";
import datasetCover from "@/assets/dataset-cover.png";
import { Image } from "@/components/ui/image";

export const Home = () => {
  return (
    <div className="w-11/12 max-w-5xl flex flex-col items-stretch gap-2">
      <div className="pt-48 h-fit">
        <Text size="t1" className="inline-block z-10">
          AI for Brain Tumor Diagnosis
        </Text>
        <div className="absolute w-screen h-124 top-0 right-0 bg-black shadow-xl">
          <img
            src={datasetCover}
            alt="Image of multiple skull MRI scans"
            className="absolute right-0 h-full object-cover"
          />
          <div className="absolute top-0 right-0 w-full h-full z-1 bg-gradient-to-l from-black to-black/30" />
        </div>
      </div>
      <Text size="h2">CS4641 Fall 2025 Team 24</Text>
      <Text size="h3">
        Ritika Calyanakoti, Aaron Fan, Hannah Hsiao, Aaron Hung, Sagarika
        Srinivasan
      </Text>
      <Text size="h2" id="overview" className="pt-8">
        Overview
      </Text>
      <div className="gap-12 clear-both">
        <Image
          float
          src={pituitaryImage}
          alt="MRI scan of human head with pituitary tumor"
          className="h-full max-h-96 rounded-lg object-contain shadow-sm"
        />
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
      </div>
    </div>
  );
};
