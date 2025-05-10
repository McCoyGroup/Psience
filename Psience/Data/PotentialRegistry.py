from McUtils.ExternalPrograms import WebAPIConnection, GitHubReleaseManager, ReleaseZIPManager
import os, importlib, zipfile, sys

__all__ = [
    "PotentialRegistryAPI",
]

class PotentialReleaseZIPManager(ReleaseZIPManager):
    default_resource_name = 'potentials'
    location_env_var = 'POTENTIAL_REGISTRY_DIR'
    use_temporary = False


class PotentialRegistryAPI(GitHubReleaseManager):
    request_base = 'https://api.github.com/'

    def __init__(self, token=None, request_delay_time=None, release_manager=None, **opts):
        if release_manager is None:
            release_manager = ReleaseZIPManager()
        self.release_manager = release_manager
        self.update_existing_releases()
        super().__init__({'header':'apikey', 'value':token}, request_delay_time=request_delay_time, **opts)

    potential_registry_org = 'Potential-Registry'
    blacklist_repos = [
        '.github'
    ]
    def list_repos(self, owner=None):
        if owner is None:
            owner = self.potential_registry_org
        return super().list_repos(owner)
    def list_potentials(self):
        return self.list_repos()

    def list_releases(self, repo_or_owner, repo=None):
        if repo is None:
            repo = repo_or_owner
            repo_or_owner = self.potential_registry_org
        return super().list_releases(repo_or_owner, repo)
    def latest_release(self,repo_or_owner, repo=None):
        if repo is None:
            repo = repo_or_owner
            repo_or_owner = self.potential_registry_org
        return super().latest_release(repo_or_owner, repo)

    def get_release_list(self, repo_or_owner, repo=None, update=None):
        if repo is None:
            repo = repo_or_owner
            repo_or_owner = self.potential_registry_org
        return super().get_release_list(repo_or_owner, repo, update=update)

    def _import_from(self, name, release_file, safety_check=True):
        if name in sys.modules:
            if safety_check:
                mod = sys.modules[name]
                loc = self.release_manager.location
                if not mod.__file__.startswith(loc):
                    raise ValueError(
                        f"module {mod} which would overwritten by loading {name} from {release_file} was not loaded from {loc}"
                    )
            del sys.modules[name]
        sys.path.insert(0, release_file)
        try:
            return importlib.import_module(name)
        finally:
            sys.path.pop(0)
    def _get_potential_module_version(self, name, version, update=None):
        releases = self.get_release_list(name, update=update)
        if version not in releases:
            raise ValueError(f"couldn't find release version {ReleaseZIPManager.make_semver(version)} for potential {name}")
        release = releases[version]
        if not os.path.exists(release):
            release = self.release_manager.get_resource(release, load_resource=False)

        return self._import_from(name, release)
    def get_potential(self, name, update=None, version=None):
        if version is None:
            releases = self.get_release_list(name, update=update)
            version = list(sorted(releases.keys()))[-1]
        if isinstance(version, str):
            version = self.release_manager.parse_semver(version)
        else:
            version = tuple(version)
        return self._get_potential_module_version(name, version, update=False)
