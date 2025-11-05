import { cn } from "@/lib/utils";
import type { ComponentProps } from "react";
import { Text } from "./text";

export interface ImageProps extends ComponentProps<"img"> {
  float?: boolean;
}

export const Image = ({
  float = false,
  className,
  alt,
  ...props
}: ImageProps) => {
  return (
    <div
      className={cn(
        `image-container flex flex-col items-center`,
        float ? "float-right" : "float-none inline-flex",
      )}
    >
      <img className={cn(className)} alt={alt} {...props} />
      <Text size="c" className="text-center max-w-11/12">
        {alt}
      </Text>
    </div>
  );
};
