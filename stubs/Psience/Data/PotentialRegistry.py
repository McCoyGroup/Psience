from McUtils.ExternalPrograms import WebAPIConnection, GitHubReleaseManager, ReleaseZIPManager
import os, importlib, zipfile, sys
__all__ = ['PotentialRegistryAPI']

class PotentialReleaseZIPManager(ReleaseZIPManager):
    default_resource_name = 'potentials'
    location_env_var = 'POTENTIAL_REGISTRY_DIR'
    use_temporary = False

class PotentialRegistryAPI(GitHubReleaseManager):
    request_base = 'https://api.github.com/'

    def __init__(self, token=None, request_delay_time=None, release_manager=None, **opts):
        ...
    potential_registry_org = 'Potential-Registry'
    blacklist_repos = ['.github']

    def list_repos(self, owner=None):
        ...

    def list_potentials(self):
        ...

    def list_releases(self, repo_or_owner, repo=None):
        ...

    def latest_release(self, repo_or_owner, repo=None):
        ...

    def get_release_list(self, repo_or_owner, repo=None, update=None):
        ...

    def _import_from(self, name, release_file, safety_check=True):
        ...

    def _get_potential_module_version(self, name, version, update=None):
        ...

    def get_potential(self, name, update=None, version=None):
        ...