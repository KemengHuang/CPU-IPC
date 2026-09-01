import argparse
from pathlib import Path

import cv2


REPOSITORY_ROOT = Path(__file__).resolve().parent.parent


def image_files(directory: Path) -> list[Path]:
    return sorted(path for path in directory.iterdir() if path.is_file())


def generate_video(image_directory: Path, output_file: Path, fps: float) -> None:
    images = image_files(image_directory)
    if not images:
        raise RuntimeError(f"no images found in {image_directory}")

    first_frame = cv2.imread(str(images[0]))
    if first_frame is None:
        raise RuntimeError(f"failed to read {images[0]}")

    height, width = first_frame.shape[:2]
    output_file.parent.mkdir(parents=True, exist_ok=True)
    writer = cv2.VideoWriter(
        str(output_file),
        cv2.VideoWriter_fourcc(*"DIVX"),
        fps,
        (width, height),
    )
    if not writer.isOpened():
        raise RuntimeError(f"failed to open video writer for {output_file}")

    try:
        for image in images:
            frame = cv2.imread(str(image))
            if frame is None:
                raise RuntimeError(f"failed to read {image}")
            if frame.shape[:2] != (height, width):
                raise RuntimeError(f"image dimensions differ: {image}")
            writer.write(frame)
    finally:
        writer.release()


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a video from viewer screenshots.")
    parser.add_argument(
        "--images",
        type=Path,
        default=REPOSITORY_ROOT / "Output" / "saveScreen",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY_ROOT / "Output" / "video.mp4",
    )
    parser.add_argument("--fps", type=float, default=10.0)
    arguments = parser.parse_args()
    generate_video(arguments.images, arguments.output, arguments.fps)


if __name__ == "__main__":
    main()
