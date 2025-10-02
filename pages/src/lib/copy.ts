import { toast } from "sonner";

export const copy = async (textToCopy: string) => {
  try {
    await navigator.clipboard.writeText(textToCopy);
    toast('URL Copied!')
  } catch (err) {
    console.error('Failed to copy text: ', err);
  }
};
