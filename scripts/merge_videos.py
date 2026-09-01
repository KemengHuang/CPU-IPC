import argparse
from pathlib import Path

import cv2


REPOSITORY_ROOT = Path(__file__).resolve().parent.parent


def image_files(directory: Path) -> list[Path]:
    return sorted(path for path in directory.iterdir() if path.is_file())


def merge_videos(
    left_directory: Path,
    right_directory: Path,
    output_file: Path,
    fps: float,
) -> None:
    left_images = image_files(left_directory)
    right_images = image_files(right_directory)
    if not left_images or len(left_images) != len(right_images):
        raise RuntimeError("the two screenshot directories must have equal non-zero counts")

    first_left = cv2.imread(str(left_images[0]))
    first_right = cv2.imread(str(right_images[0]))
    if first_left is None or first_right is None:
        raise RuntimeError("failed to read the first screenshot pair")
    if first_left.shape != first_right.shape:
        raise RuntimeError("the two screenshot streams must have identical dimensions")

    height, single_width = first_left.shape[:2]
    output_file.parent.mkdir(parents=True, exist_ok=True)
    writer = cv2.VideoWriter(
        str(output_file),
        cv2.VideoWriter_fourcc(*"DIVX"),
        fps,
        (2 * single_width, height),
    )
    if not writer.isOpened():
        raise RuntimeError(f"failed to open video writer for {output_file}")

    try:
        for left_path, right_path in zip(left_images, right_images):
            left = cv2.imread(str(left_path))
            right = cv2.imread(str(right_path))
            if left is None or right is None or left.shape != right.shape:
                raise RuntimeError(f"invalid screenshot pair: {left_path}, {right_path}")
            writer.write(cv2.hconcat([left, right]))
    finally:
        writer.release()


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a side-by-side comparison video.")
    parser.add_argument(
        "--left",
        type=Path,
        default=REPOSITORY_ROOT / "Output" / "saveScreen1",
    )
    parser.add_argument(
        "--right",
        type=Path,
        default=REPOSITORY_ROOT / "Output" / "saveScreen2",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY_ROOT / "Output" / "comparison.mp4",
    )
    parser.add_argument("--fps", type=float, default=30.0)
    arguments = parser.parse_args()
    merge_videos(
        arguments.left,
        arguments.right,
        arguments.output,
        arguments.fps,
    )


if __name__ == "__main__":
    main()
