from McUtils.ExternalPrograms import WebAPIConnection, WebResourceManager
import os, importlib, zipfile, sys

__all__ = [
    "PotentialRegistryAPI",
]

class ReleaseZIPManager(WebResourceManager):
    default_resource_name = 'potentials'
    location_env_var = 'POTENTIAL_REGISTRY_DIR'
    use_temporary = False

    @classmethod
    def parse_semver(cls, version_string):
        return tuple(
            int(t)
            for t in version_string.rsplit('v', 1)[-1].split('.')
        )

    @classmethod
    def make_semver(cls, version):
        return 'v' + '.'.join(str(v) for v in version)

    @classmethod
    def parse_name_version(cls, filename):
        name, version = filename.rsplit('-', 1)[0].rsplit('v', 1)
        return name, cls.parse_semver(version)

    def list_resources(self):
        base_list = super().list_resources()
        resource_list = {}
        for filename, filepath in base_list.items():
            try:
                name, version = self.parse_name_version(filename)
            except ValueError:
                pass
            else:
                if name not in resource_list: resource_list[name] = {}
                resource_list[name][version] = filepath
        return resource_list
    def save_resource(self, loc, val):
        super().save_resource(loc, val)
        with zipfile.ZipFile(loc) as zip_extractor:
            dirs = zip_extractor.namelist()
            github_dirname = dirs[0]
            cwd = os.getcwd()
            try:
                os.chdir(self.location)
                zip_extractor.extractall(members=[d for d in dirs if d.startswith(github_dirname)])
                try:
                    os.remove(loc)
                except:
                    raise
                else:
                    os.rename(github_dirname, loc)
            finally:
                os.chdir(cwd)


class PotentialRegistryAPI(WebAPIConnection):
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
    def list_potentials(self):
        return {
            r['name']:r
            for r in self.get('orgs', self.potential_registry_org, 'repos')
            if r['name'] not in self.blacklist_repos
        }

    def list_releases(self, repo):
        return {
            r['tag_name']:r
            for r in self.get('repos', self.potential_registry_org, repo, 'releases')
        }
    def latest_release(self, repo):
        return self.get('repos', self.potential_registry_org, repo, 'releases', 'latest')

    # release_loc = None
    # def download_release(self, release_dict, where=None):
    #     with urllib.request.urlopen(release_dict['zipball_url']) as f:
    #         zip = f.read()

    release_cache = {}
    def update_existing_releases(self):
        release_list = self.release_manager.list_resources()
        for name, version_list in release_list.items():
            if name not in self.release_cache: self.release_cache[name] = {}
            self.release_cache[name].update(version_list)

    def get_release_list(self, name, update=None):
        if update or (update is None and name not in self.release_cache):
            self.release_cache[name] = {
                ReleaseZIPManager.parse_semver(k): v['zipball_url']
                for k, v in self.list_releases(name).items()
            }
        elif name not in self.release_cache:
            self.update_existing_releases()
        if name not in self.release_cache:
            raise ValueError(f"couldn't find any potential releases for {name}")
        return self.release_cache[name]

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
            version = ReleaseZIPManager.parse_semver(version)
        else:
            version = tuple(version)
        return self._get_potential_module_version(name, version, update=False)
