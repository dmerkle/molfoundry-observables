"""Execute a manifest's declared counting adapter on explicit native inputs."""
from dataclasses import dataclass
from mfobs.api import ModelInputs
from mfobs.inputs import ModelManifest


@dataclass(frozen=True)
class ManifestBackend:
    manifest: ModelManifest

    def run(self, inputs: ModelInputs, **options) -> dict:
        from mfobs.core.worklist import run
        inputs.require(self.manifest.specification())
        return run(self.manifest, **options)
