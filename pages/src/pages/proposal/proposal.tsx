import { Divider } from "@/components/ui/divider";
import { Text } from "@/components/ui/text";

export const Proposal = () => {
  return (
    <div className="flex flex-col w-11/12 max-w-5xl items-stretch">
      <Text size="t2" className="mb-8">
        Proposal
      </Text>
      <Text size="h2" id="introduction">
        Introduction
      </Text>
      <Divider />
      <Text>
        Brain tumors are among the most aggressive diseases, requiring early and
        accurate diagnosis for effective treatment. Analysis of brain MRI scans
        is the healthcare standard for brain tumor diagnosis and simultaneously
        distinguishes between various types of cancers and limits patient
        exposure to radiation. (CITATION). MRIs also measure treatment response.
        It is useful for oncologists to decide when to stop treatments, if other
        methods should be incorporated, etc. (CITATION). That being said, MRI
        interpretations by clinicians are not fully accurate. The paper “The ins
        and outs of errors in oncology imaging” discuss radiologist errors in
        interpreting cancer imaging such as misreading scans or having their own
        cognitive biases. Such scenarios motivate healthcare professionals to
        seek forms of diagnosis that further accuracy (CITATION). The complexity
        of brain tumors further complicates diagnosis. The tumor imaging space
        is multi-dimentional with clinicians relying on MRI to differentiate
        among numerous tumor classifications. Examples include Gliomas,
        Meningiomas, Lymphomas, Pineal tumors, etc. These cancers spread to
        different parts of the brain with some affecting the cerebrospinal
        fluid, the soft tissue, etc. and others affecting the early and mid
        brain structures (CITATION). With so many tumor types and affected
        regions, it is key to have a reference tumor aggressiveness/intensity
        scale. The World Health Organization's grading system has been the
        widely accepted standard which classifies the central nervous system
        tumors from Grade I (least malignant) to Grade IV (most malignant)
        (CITATION). While the WHO grading system is a standardized framework,
        its application still relies heavily on radiologist interpretation of
        MRI scans. This dependence introduces potential error, creating a need
        for computational methods that can assist in improving diagnostic
        consistency and accuracy. In this context, artificial intelligence (AI)
        and machine learning (ML) approaches have emerged as promising tools in
        medical imaging. AI/ML applications in healthcare are now increasingly
        common with detecting rates of skin cancer, drug development
        AI-Monitoring sensors (CITATION). This project provides an exciting
        opportunity to enhance diagnostic accuracy and overall clinical
        decision-making.
      </Text>
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
      <Text size="h2" id="references">
        Methods
      </Text>
      <Divider />
      <Text size="h2" id="references">
        Results and Discussion
      </Text>
      <Divider />
      <Text></Text>
      <Text size="h2" id="references">
        References
      </Text>
      <Divider />
      <ol>
        <li>
          <Text>
            A. Iannessi, H. Beaumont, C. Aguillera, F. Nicol, and A.-S.
            Bertrand, “The ins and outs of errors in oncology imaging: the DAC
            framework for radiologists,” Frontiers in Oncology, vol. 14, Oct.
            2024, doi: https://doi.org/10.3389/fonc.2024.1402838.
          </Text>
        </li>
        <li>
          <Text>
            D. N. Louis et al., “The 2021 WHO Classification of Tumors of the
            Central Nervous System: a summary,” Neuro-Oncology, vol. 23, no. 8,
            Jun. 2021, doi: https://doi.org/10.1093/neuonc/noab106.
          </Text>
        </li>
        <li>
          <Text>
            M. Martucci et al., “Magnetic Resonance Imaging of Primary Adult
            Brain Tumors: State of the Art and Future Perspectives,”
            Biomedicines, vol. 11, no. 2, p. 364, Feb. 2023, doi:
            https://doi.org/10.3390/biomedicines11020364.
          </Text>
        </li>
        <li>
          <Text>
            N. Melarkode, K. Srinivasan, S. M. Qaisar, and P. Plawiak,
            “AI-Powered Diagnosis of Skin Cancer: A Contemporary Review, Open
            Challenges and Future Research Directions,” Cancers, vol. 15, no. 4,
            p. 1183, Jan. 2023, doi: https://doi.org/10.3390/cancers15041183.
          </Text>
        </li>
        <li>
          <Text>
            S. R. Chandana, S. Movva, M. Arora, and T. Singh, “Primary Brain
            Tumors in Adults,” American Family Physician, vol. 77, no. 10, pp.
            1423–1430, May 2008, Available:
            https://www.aafp.org/pubs/afp/issues/2008/0515/p1423.html
          </Text>
        </li>
      </ol>
    </div>
  );
};
