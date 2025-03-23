"""Utilities for working with lists."""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path
import tempfile
from typing import Optional, Sequence
import uuid

from ase.atoms import Atoms
from ase.io import write


def render_atoms_trajectory(
    trajectory: Sequence[Atoms],
    output_dir: str,
    image_config: Optional[dict[str, str]] | None = None,
    video_config: Optional[dict[str, str]] | None = None,
) -> None:
    """Generate a video of the atomic trajectory using ASE and state final image

    Parameters
    ----------
    trajectory: Sequence of ASE Atoms objects representing trajectory frames
    output_dir: Directory to save the video/image files
    image_config
        - rotation: Rotation string (e.g. '10z,-80x')
        - scale: Scale factor for image (default: 20)
        - show_unit_cell: Show unit cell (default: False)
        - format: Image format (default: 'png')
    video_config
        - rotation: Rotation string (e.g. '10z,-80x')
        - scale: Scale factor for image (default: 20)
        - fps: Frames per second in output video (default: 5)
        - width: Width of output video (default: 800)
        - height: Height of output video (default: 600)

    Returns
    -------
    Void
    """
    if output_dir is None:
        raise ValueError("output_dir must be specified")

    if len(trajectory) < 2:
        raise ValueError("Trajectory must contain at multiple frames")

    # Set/get defaults
    image_config = {
        "style": "cartoon",
        "colorscheme": "bwr",
        "ray": True,
        "width": 1024,
        "height": 768,
    } | (image_config or {})
    video_config = {"fps": 5, "width": 800, "height": 600} | (video_config or {})

    # Extract and validate config
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(exist_ok=True, parents=True)

    # Generate unique filename based on trajectory properties
    first_atoms = trajectory[0]
    chemical_formula = first_atoms.get_chemical_formula(mode="hill")
    num_atoms = len(first_atoms)
    num_frames = len(trajectory)
    structure_hash = hashlib.md5(
        str(first_atoms.positions + trajectory[-1].positions).encode()
    ).hexdigest()[:8]
    unique_id = str(uuid.uuid4())[:8]
    output_filename_prefix = (
        f"{chemical_formula}_{num_atoms}atoms_{structure_hash}_{unique_id}"
    )

    # VIDEO: TRAJECTORY
    width = video_config.get("width", 800)
    height = video_config.get("height", 600)
    width += width % 2
    height += height % 2

    # Use a temporary directory for frames
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = Path(temp_dir)

        # Render frames in temporary directory
        for i, atoms in enumerate(trajectory):
            frame_path = temp_dir_path / f"frame_{i:04d}.png"
            write(
                str(frame_path),
                atoms,
                format="png",
                show_unit_cell=0,
                rotation=video_config.get("rotation", "0z"),
                scale=video_config.get("scale", 20),
            )

        # Generate video from frames
        output_video_path = (
            output_dir / f"{output_filename_prefix}_{num_frames}frames.mp4"
        )
        cmd_str = (
            f"ffmpeg -y -framerate {video_config.get('fps', 5)} "
            f"-i {temp_dir_path}/frame_%04d.png "
            f"-vf scale={width}:{height} "
            f"-c:v libx264 -pix_fmt yuv420p {output_video_path}"
        )
        subprocess.run(cmd_str, shell=True, check=True)

    # IMAGE: FINAL
    output_final_frame_path = output_dir / f"{output_filename_prefix}_final.png"
    write(
        str(output_final_frame_path),
        trajectory[-1],
        format=image_config.get("format", "png"),
        show_unit_cell=image_config.get("show_unit_cell", 0),
        rotation=image_config.get("rotation", "0z"),
        scale=image_config.get("scale", 20),
    )
