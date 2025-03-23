"""Utilities for working with lists."""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path
import shutil
from typing import Sequence
import uuid

from ase.atoms import Atoms
from ase.io import write


def render_atoms(
    trajectory: Sequence[Atoms],
    renders_dir: str,
    render_images: bool = False,
    render_video: bool = False,
    image_config: dict[str, str] | None = None,
    video_config: dict[str, str] | None = None,
) -> None:
    """Generate a video of the atomic trajectory using ASE.

    Parameters
    ----------
    trajectory
        Sequence of ASE Atoms objects representing the trajectory frames
    renders_dir: Directory to save the video file
    filename_suffix: Suffix for the output file name
    image_config
        Configuration options for image generation, including:
        - rotation: Rotation string (e.g. '10z,-80x')
        - scale: Scale factor for image (default: 20)
        - show_unit_cell: Show unit cell (default: False)
        - format: Image format (default: 'png')
    video_config
        Configuration options for video generation, including:
        - rotation: Rotation string (e.g. '10z,-80x')
        - scale: Scale factor for image (default: 20)
        - fps: Frames per second in output video (default: 5)
        - width: Width of output video (default: 800)
        - height: Height of output video (default: 600)

    Returns
    -------
    Path
        Path to the generated video file
    """
    if render_images != True and render_video != True:
        return

    # Extract and validate config
    output_dir = Path(renders_dir).resolve()
    output_dir.mkdir(exist_ok=True, parents=True)

    # Create unique frames directory
    unique_id = str(uuid.uuid4())[:8]  # Use first 8 chars of UUID
    frames_dir = output_dir / f"frames_{unique_id}"
    frames_dir.mkdir(exist_ok=True, parents=True)

    # Generate unique filename based on trajectory properties
    first_atoms = trajectory[0]
    chemical_formula = first_atoms.get_chemical_formula(mode="hill")
    num_atoms = len(first_atoms)
    num_frames = len(trajectory)
    structure_hash = hashlib.md5(
        str(first_atoms.positions + trajectory[-1].positions).encode()
    ).hexdigest()[:8]

    # FRAMES: INITIAL (if more than 1 frame)
    if render_images == True and len(trajectory) > 1:
        output_initial_frame_path = (
            output_dir
            / f"{chemical_formula}_{num_atoms}atoms_{structure_hash}_{unique_id}_initial.png"
        )
        write(
            str(output_initial_frame_path),
            trajectory[0],
            format=image_config.get("format", "png"),
            show_unit_cell=image_config.get("show_unit_cell", 0),
            rotation=image_config.get("rotation", "0z"),
            scale=image_config.get("scale", 20),
        )

    # FRAMES: TRAJECTORY (if more than 1 frame)
    if render_video == True and len(trajectory) > 1:
        # Set default dimensions that are guaranteed to be even
        width = video_config.get("width", 800)  # Default width
        height = video_config.get("height", 600)  # Default height
        width += width % 2
        height += height % 2
        # Generate frames
        for i, atoms in enumerate(trajectory):
            frame_path = frames_dir / f"frame_{i:04d}.png"
            # ASE's write function can directly output PNG files
            write(
                str(frame_path),
                atoms,
                format="png",
                show_unit_cell=0,
                rotation=video_config.get("rotation", "0z"),
                scale=video_config.get("scale", 20),
            )
            if not frame_path.exists():
                raise RuntimeError(f"Failed to save frame {i}")

        # Render File
        output_file = (
            f"{chemical_formula}_{num_atoms}atoms_"
            f"{structure_hash}_"
            f"{unique_id}_"
            f"{num_frames}frames.mp4"
        )
        output_path = output_dir / output_file
        cmd_str = (
            f"ffmpeg -y -framerate {video_config.get('fps', 5)} "
            f"-i {frames_dir}/frame_%04d.png "
            f"-vf scale={width}:{height} "  # Force even dimensions
            f"-c:v libx264 -pix_fmt yuv420p {output_path}"
        )
        subprocess.run(cmd_str, shell=True, check=True)

        # Clean up frames directory after video is created
        shutil.rmtree(frames_dir, ignore_errors=True)

    # FRAMES: FINAL (always)
    output_final_frame_path = (
        output_dir
        / f"{chemical_formula}_{num_atoms}atoms_{structure_hash}_{unique_id}_final.png"
    )
    write(
        str(output_final_frame_path),
        trajectory[-1],
        format=image_config.get("format", "png"),
        show_unit_cell=image_config.get("show_unit_cell", 0),
        rotation=image_config.get("rotation", "0z"),
        scale=image_config.get("scale", 20),
    )
