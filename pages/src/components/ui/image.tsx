import { cn } from "@/lib/utils";
import type { ComponentProps } from "react";
import { Text } from "./text";

export const Image = ({ className, alt, ...props }: ComponentProps<"img">) => {
  return (
    <div className="float-right flex flex-col items-stretch">
      <img className={cn(className)} alt={alt} {...props} />
      <Text size="c" className="text-center">
        {alt}
      </Text>
    </div>
  );
};
