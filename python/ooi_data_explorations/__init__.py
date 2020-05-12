import pkgutil
import yaml

from munch import Munch

m2m_urls = pkgutil.get_data(__name__, 'm2m_urls.yml')
M2M_URLS = Munch.fromDict(yaml.safe_load(m2m_urls))
