import datetime
import os
import time
from typing import Callable

import requests
from dotenv import load_dotenv


def telegram_logger(chat_id: int) -> Callable:
    def decorator(func):
        load_dotenv("secure.env")
        tg_url = "https://api.telegram.org/bot"

        def wrapper(*args, **kwargs):
            start = time.time()
            try:
                func(*args, **kwargs)
            except Exception as error:
                message = f"Function `{func.__name__}` failed with an exception:\n\n`{type(error).__name__}: {error}`"
            else:
                end = time.time()
                execution_time = end - start

                time_obj = datetime.datetime.fromtimestamp(execution_time, datetime.UTC)
                if execution_time // 86400 < 1:
                    formatted_time = time_obj.strftime("%H:%M:%S.%f")
                else:
                    formatted_time = time_obj.strftime("%-j days, %H:%M:%S")

                message = f"Function `{func.__name__}` finished in `{formatted_time}`"

            requests.post(tg_url + os.getenv("TG_API_TOKEN") + "/sendMessage",
                          params={"chat_id": chat_id, "parse_mode": "MarkdownV2", "text": message})

        return wrapper

    return decorator
